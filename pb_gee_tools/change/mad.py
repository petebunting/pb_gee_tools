import numpy
import ee


def _covarw(
    image: ee.Image,
    weights: ee.Image = None,
    scale: float = 20,
    max_pixels: float = 1e10,
):
    """Return the centered image and its weighted covariance matrix."""
    try:
        geometry = image.geometry()
        band_names = image.bandNames()
        N = band_names.length()
        if weights is None:
            weights = image.constant(1)
        weights_image = image.multiply(ee.Image.constant(0)).add(weights)
        means = (
            image.addBands(weights_image)
            .reduceRegion(
                ee.Reducer.mean().repeat(N).splitWeights(),
                scale=scale,
                maxPixels=max_pixels,
            )
            .toArray()
            .project([1])
        )
        centered = image.toArray().subtract(means)
        B1 = centered.bandNames().get(0)
        b1 = weights.bandNames().get(0)
        n_pixels = ee.Number(
            centered.reduceRegion(
                ee.Reducer.count(), scale=scale, maxPixels=max_pixels
            ).get(B1)
        )
        sum_weights = ee.Number(
            weights.reduceRegion(
                ee.Reducer.sum(), geometry=geometry, scale=scale, maxPixels=max_pixels
            ).get(b1)
        )
        covw = (
            centered.multiply(weights.sqrt())
            .toArray()
            .reduceRegion(
                ee.Reducer.centeredCovariance(),
                geometry=geometry,
                scale=scale,
                maxPixels=max_pixels,
            )
            .get("array")
        )
        covw = ee.Array(covw).multiply(n_pixels).divide(sum_weights)
        return (centered.arrayFlatten([band_names]), covw)
    except Exception as e:
        print("Error: %s" % e)
        raise e


"""
def _trans_cov_to_corr(cov):
    # Transfrom covariance matrix to correlation matrix.

    def _trunc(values, dec=3):
        # Truncate a 1-D array to dec decimal places.
        return numpy.trunc(values * 10 ** dec) / (10 ** dec)
    
    # Diagonal matrix of inverse sigmas.
    s_inv = cov.matrixDiagonal().sqrt().matrixToDiag().matrixInverse()
    # Transform.
    corr = s_inv.matrixMultiply(cov).matrixMultiply(s_inv).getInfo()
    # Truncate.
    return [list(map(_trunc, corr[i])) for i in range(len(corr))]
"""


def _geneiv(C, B):
    """Return the eignvalues and eigenvectors of the generalized eigenproblem
    C*X = lambda*B*X"""
    try:
        C = ee.Array(C)
        B = ee.Array(B)
        #  Li = choldc(B)^-1
        Li = ee.Array(B.matrixCholeskyDecomposition().get("L")).matrixInverse()
        # Solve symmetric, ordinary eigenproblem Li*C*Li^T*x = lambda*x
        Xa = Li.matrixMultiply(C).matrixMultiply(Li.matrixTranspose()).eigen()
        # Eigenvalues as a row vector.
        lambdas = Xa.slice(1, 0, 1).matrixTranspose()
        # Eigenvectors as columns.
        X = Xa.slice(1, 1).matrixTranspose()
        # Generalized eigenvectors as columns, Li^T*X
        eigenvecs = Li.matrixTranspose().matrixMultiply(X)
        return (lambdas, eigenvecs)
    except Exception as e:
        print("Error: %s" % e)
        raise e


def mad_run(
    image1: ee.Image, image2: ee.Image, scale: float = 20
) -> [ee.Image, ee.Image, ee.Image, ee.Image]:
    """The MAD transformation of two multiband images."""
    try:
        image = image1.addBands(image2)
        n_bands = image.bandNames().length().divide(2)
        centered_image, covarArray = _covarw(image, scale=scale)
        b_names = centered_image.bandNames()
        b_names1 = b_names.slice(0, n_bands)
        b_names2 = b_names.slice(n_bands)
        centered_image1 = centered_image.select(b_names1)
        centered_image2 = centered_image.select(b_names2)
        s11 = covarArray.slice(0, 0, n_bands).slice(1, 0, n_bands)
        s22 = covarArray.slice(0, n_bands).slice(1, n_bands)
        s12 = covarArray.slice(0, 0, n_bands).slice(1, n_bands)
        s21 = covarArray.slice(0, n_bands).slice(1, 0, n_bands)
        c1 = s12.matrixMultiply(s22.matrixInverse()).matrixMultiply(s21)
        b1 = s11
        c2 = s21.matrixMultiply(s11.matrixInverse()).matrixMultiply(s12)
        b2 = s22

        # Solution of generalized eigenproblems.
        lambdas, A = _geneiv(c1, b1)
        _, B = _geneiv(c2, b2)
        rhos = lambdas.sqrt().project(ee.List([1]))

        # MAD variances.
        sigma2s = rhos.subtract(1).multiply(-2).toList()
        sigma2s = ee.Image.constant(sigma2s)

        # Ensure sum of positive correlations between X and U is positive.
        tmp = s11.matrixDiagonal().sqrt()
        ones = tmp.multiply(0).add(1)
        tmp = ones.divide(tmp).matrixToDiag()
        s = (
            tmp.matrixMultiply(s11)
            .matrixMultiply(A)
            .reduce(ee.Reducer.sum(), [0])
            .transpose()
        )
        A = A.matrixMultiply(s.divide(s.abs()).matrixToDiag())

        # Ensure positive correlation.
        tmp = A.transpose().matrixMultiply(s12).matrixMultiply(B).matrixDiagonal()
        tmp = tmp.divide(tmp.abs()).matrixToDiag()
        B = B.matrixMultiply(tmp)

        # Canonical and MAD variates as images.
        centered_image1_array = centered_image1.toArray().toArray(1)
        centered_image2_array = centered_image2.toArray().toArray(1)
        U = (
            ee.Image(A.transpose())
            .matrixMultiply(centered_image1_array)
            .arrayProject([0])
            .arrayFlatten([b_names2])
        )
        V = (
            ee.Image(B.transpose())
            .matrixMultiply(centered_image2_array)
            .arrayProject([0])
            .arrayFlatten([b_names2])
        )
        MAD = U.subtract(V)

        # Chi-square image.
        Z = MAD.pow(2).divide(sigma2s).reduce(ee.Reducer.sum())
        return (U, V, MAD, Z)
    except Exception as e:
        print("Error: %s" % e)
        raise e


def _imad(current, prev):
    """Iterator function for iMAD."""

    def _imad1(current, prev):
        """Iteratively re-weighted MAD."""

        def _chi2cdf(Z, df):
            """Chi-square cumulative distribution function with df degrees of freedom."""
            return ee.Image(Z.divide(2)).gammainc(ee.Number(df).divide(2))

        image = ee.Image(ee.Dictionary(prev).get("image"))
        Z = ee.Image(ee.Dictionary(prev).get("Z"))
        allrhos = ee.List(ee.Dictionary(prev).get("allrhos"))
        n_bands = image.bandNames().length().divide(2)
        weights = _chi2cdf(Z, n_bands).subtract(1).multiply(-1)
        scale = ee.Dictionary(prev).getNumber("scale")
        niter = ee.Dictionary(prev).getNumber("niter")

        # Weighted stacked image and weighted covariance matrix.
        centeredImage, covarArray = _covarw(image, weights, scale)
        b_names = centeredImage.bandNames()
        b_names1 = b_names.slice(0, n_bands)
        b_names2 = b_names.slice(n_bands)
        centered_image1 = centeredImage.select(b_names1)
        centered_image2 = centeredImage.select(b_names2)
        s11 = covarArray.slice(0, 0, n_bands).slice(1, 0, n_bands)
        s22 = covarArray.slice(0, n_bands).slice(1, n_bands)
        s12 = covarArray.slice(0, 0, n_bands).slice(1, n_bands)
        s21 = covarArray.slice(0, n_bands).slice(1, 0, n_bands)
        c1 = s12.matrixMultiply(s22.matrixInverse()).matrixMultiply(s21)
        b1 = s11
        c2 = s21.matrixMultiply(s11.matrixInverse()).matrixMultiply(s12)
        b2 = s22

        # Solution of generalized eigenproblems.
        lambdas, A = _geneiv(c1, b1)
        _, B = _geneiv(c2, b2)
        rhos = lambdas.sqrt().project(ee.List([1]))
        # Test for convergence.
        lastrhos = ee.Array(allrhos.get(-1))
        done = (
            rhos.subtract(lastrhos)
            .abs()
            .reduce(ee.Reducer.max(), ee.List([0]))
            .lt(ee.Number(0.0001))
            .toList()
            .get(0)
        )
        allrhos = allrhos.cat([rhos.toList()])

        # MAD variances.
        sigma2s = rhos.subtract(1).multiply(-2).toList()
        sigma2s = ee.Image.constant(sigma2s)
        # Ensure sum of positive correlations between X and U is positive.
        tmp = s11.matrixDiagonal().sqrt()
        ones = tmp.multiply(0).add(1)
        tmp = ones.divide(tmp).matrixToDiag()
        s = (
            tmp.matrixMultiply(s11)
            .matrixMultiply(A)
            .reduce(ee.Reducer.sum(), [0])
            .transpose()
        )
        A = A.matrixMultiply(s.divide(s.abs()).matrixToDiag())

        # Ensure positive correlation.
        tmp = A.transpose().matrixMultiply(s12).matrixMultiply(B).matrixDiagonal()
        tmp = tmp.divide(tmp.abs()).matrixToDiag()
        B = B.matrixMultiply(tmp)

        # Canonical and MAD variates.
        centered_image1_array = centered_image1.toArray().toArray(1)
        centered_image2_array = centered_image2.toArray().toArray(1)
        U = (
            ee.Image(A.transpose())
            .matrixMultiply(centered_image1_array)
            .arrayProject([0])
            .arrayFlatten([b_names1])
        )
        V = (
            ee.Image(B.transpose())
            .matrixMultiply(centered_image2_array)
            .arrayProject([0])
            .arrayFlatten([b_names2])
        )
        iMAD = U.subtract(V)

        # Chi-square image.
        Z = iMAD.pow(2).divide(sigma2s).reduce(ee.Reducer.sum())
        return ee.Dictionary(
            {
                "done": done,
                "scale": scale,
                "niter": niter.add(1),
                "image": image,
                "allrhos": allrhos,
                "Z": Z,
                "iMAD": iMAD,
            }
        )

    done = ee.Number(ee.Dictionary(prev).get("done"))
    return ee.Algorithms.If(done, prev, _imad1(current, prev))


def run_imad(
    aoi: ee.Geometry,
    image1: ee.Image,
    image2: ee.Image,
    scale: float = 20,
    max_iter: int = 100,
) -> [ee.Image, ee.String, ee.Image, ee.Number]:
    try:
        N = image1.bandNames().length().getInfo()
        imad_names = ["iMAD" + str(i + 1) for i in range(N)]
        imad_names.append("Z")
        # Maximum iterations.
        input_list = ee.List.sequence(1, max_iter)
        first = ee.Dictionary(
            {
                "done": 0,
                "scale": scale,
                "niter": ee.Number(0),
                "image": image1.addBands(image2),
                "allrhos": [ee.List.sequence(1, N)],
                "Z": ee.Image.constant(0),
                "iMAD": ee.Image.constant(0),
            }
        )
        # Iteration.
        result = ee.Dictionary(input_list.iterate(_imad, first))
        # Retrieve results.
        iMAD = ee.Image(result.get("iMAD")).clip(aoi)
        rhos = ee.String.encodeJSON(ee.List(result.get("allrhos")).get(-1))
        Z = ee.Image(result.get("Z"))
        niter = result.getNumber("niter")
        # create iMAD and Z as a singe image, including rhos and number of iterations in properties.
        # imad_img = (
        #    ee.Image.cat(iMAD, Z).rename(imad_names).set("rhos", rhos, "niter", niter)
        # )
        return (iMAD, rhos, Z, niter)
    except Exception as e:
        print("Error: %s" % e)
        raise e
