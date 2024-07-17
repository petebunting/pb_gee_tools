#!/usr/bin/env python

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
# Purpose:  Get access to datasets.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 17/07/2025
# Version: 1.0
#
# History:
# Version 1.0 - Created.

import datetime


def _do_dates_overlap(
    s1_date: datetime.datetime,
    e1_date: datetime.datetime,
    s2_date: datetime.datetime,
    e2_date: datetime.datetime,
):
    """
    A function which checks if two dated periods overlap.

    :param s1_date: start of period 1
    :param e1_date: end of period 1
    :param s2_date: start of period 2
    :param e2_date: end of period 2

    :return: boolean - True the two periods overlap, False otherwise
    """
    latest_start = max(s1_date, s2_date)
    earliest_end = min(e1_date, e2_date)
    delta = (earliest_end - latest_start).days + 1
    overlap = max(0, delta)
    return overlap > 0
