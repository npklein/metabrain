"""
File:         df_utilities.py
Created:      2020/03/12
Last Changed: 2020/03/13
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

# Local application imports.
from src.utilities import get_basename


def load_dataframe(inpath, header, index_col, sep="\t", low_memory=True):
    df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                     low_memory=low_memory)
    print("\tLoaded dataframe: {} with shape: {}".format(get_basename(inpath),
                                                         df.shape))
    return df


def save_dataframe(df, outpath, header, index, sep="\t"):
    compression = 'infer'
    if outpath.endswith('.gz'):
        compression = 'gzip'

    df.to_csv(outpath, sep=sep, index=index, header=header,
              compression=compression)
    print("\tSaved dataframe: {} with shape: {}".format(get_basename(outpath),
                                                        df.shape))
