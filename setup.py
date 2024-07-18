#!/usr/bin/env python
"""
Setup script for pb_gee_tools. Use like this for Unix:

$ pip install .

"""
# This file is part of 'pb_gee_tools'
#
# Copyright 2024 Pete Bunting
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#
# Purpose:  install software.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 17/07/2025
# Version: 1.0
#
# History:
# Version 1.0 - Created.

import glob

import setuptools

import pb_gee_tools

setuptools.setup(
    name="pb_gee_tools",
    version=pb_gee_tools.PB_GEE_TOOLS_VERSION,
    description="A module with tools for supporting the use of the Google Earth Engine Python API and reducing work for repetitive or common tasks.",
    author="Pete Bunting",
    author_email="petebunting@mac.com",
    include_package_data=True,
    #scripts=glob.glob("bin/*.py"),
    packages=["pb_gee_tools", "pb_gee_tools/change"],
    license="LICENSE.txt",
    install_requires=["earthengine-api"],
    url="https://github.com/petebunting/pb_gee_tools",
    classifiers=[
        "Intended Audience :: Developers",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
)
