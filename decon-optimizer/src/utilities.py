"""
File:         utilities.py
Created:      2020/12/21
Last Changed:
Author(s):    M.Vochteloo

Copyright (C) 2020 M.Vochteloo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License can be found in the LICENSE file in the
root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
"""

# Standard imports.

# Third party imports.
import pandas as pd
from scipy import stats

# Local application imports.


def force_normal(df, value_axis=1, log=None):
    work_df = df.copy()

    if value_axis == 0:
        work_df = work_df.T
    elif value_axis == 1:
        pass
    else:
        if log is not None:
            log.error("Unexpected axis in force normal function.")
        exit()

    data = []
    for index, row in work_df.iterrows():
        data.append(force_normal_series(row))

    normal_df = pd.DataFrame(data, index=work_df.index, columns=work_df.columns)

    if value_axis == 0:
        normal_df = normal_df.T

    return normal_df


def force_normal_series(s, as_series=False):
    normal_s = stats.norm.ppf((s.rank(ascending=True) - 0.5) / s.size)
    if as_series:
        return pd.Series(normal_s, index=s.index)
    else:
        return normal_s