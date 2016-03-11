#!/usr/bin/python
"""ncsm_single_calc.py

To run as a script:

    $ ncsm_single_calc.py [-[Ff]] Z A [Aeff [nhw [n1 n2] | n1 n2]]

In the current directory, run NCSD for a single element with a single A
value and Aeff value.

If arguments are preceded by -f or -F, force recalculation of NCSD, even if
    outfiles already exist. Otherwise, calculations are not repeated.
If 2 arguments are given, these are Z A. Aeff is taken to be equal to A.
If 3 arguments are given, these are Z A Aeff.
If 4 arguments are given, these are Z A Aeff nhw.
If 5 arguments are given, these are Z A Aeff n1 n2.
If 6 arguments are given, these are Z A Aeff nhw n1 n2.
"""
from __future__ import division
from os import path
from sys import argv
from InvalidNumberOfArgumentsException import InvalidNumberOfArgumentsException
from ncsm_vce_calc import get_name, make_base_directories, make_mfdp_files
from ncsm_vce_calc import truncate_spaces, do_ncsd
from ncsm_vce_calc import NHW, N1, N2
from ncsm_vce_calc import PATH_RESULTS, PATH_TEMPLATES, DIR_FMT_NUC
from ncsm_vce_calc import REGEX_TBME, FNAME_FMT_NCSD_OUT, FNAME_MFDP


def ncsm_single_calculation(
        z, a, aeff,
        nhw=NHW, n1=N1, n2=N2,
        _path_results=PATH_RESULTS,
        _path_temp=PATH_TEMPLATES,
        _dir_fmt_nuc=DIR_FMT_NUC,
        _fname_regex_tbme=REGEX_TBME,
        _fname_fmt_ncsd_out=FNAME_FMT_NCSD_OUT,
        _fname_mfdp=FNAME_MFDP,
        force=False):
    # get directory path
    dirpath_nuc = path.join(
        _path_results,
        _dir_fmt_nuc % (get_name(z=z), a, aeff, nhw, n1, n2))
    outfile_ncsd = path.join(
        dirpath_nuc,
        _fname_fmt_ncsd_out % (get_name(z=z), a, aeff, nhw, n1, n2))
    # make directory
    make_base_directories(a_values=[a], presc=[aeff],
                          results_path=_path_results,
                          a_dirpaths_map={a: dirpath_nuc})

    # ncsm calculations: make mfdp file, perform truncation, and run NCSD
    make_mfdp_files(z=z, a_range=[a], a_presc=[aeff],
                    a_dirpath_map={a: dirpath_nuc},
                    a_outfile_map={a: outfile_ncsd},
                    n_hw=nhw, n_1=n1, n_2=n2,
                    mfdp_name=_fname_mfdp)
    truncate_spaces(n1=n1, n2=n2, dirpaths=[dirpath_nuc],
                    path_temp=_path_temp,
                    tbme_name_regex=_fname_regex_tbme)
    do_ncsd(a_values=[a], presc=[aeff],
            a_dirpaths_map={a: dirpath_nuc},
            a_outfile_map={a: outfile_ncsd},
            force=force)


if __name__ == "__main__":
    force = False
    if '-' in argv[1]:
        if 'f' in argv[1] or 'F' in argv[1]:
            force = True
        user_args = argv[2:]
    else:
        user_args = argv[1:]
    if len(user_args) == 2:
        z0, a0 = [int(x) for x in user_args]
        ncsm_single_calculation(z=z0, a=a0, aeff=a0)
    elif len(user_args) == 3:
        z0, a0, aeff0 = [int(x) for x in user_args]
        ncsm_single_calculation(z=z0, a=a0, aeff=aeff0)
    elif len(user_args) == 4:
        z0, a0, aeff0, nhw0 = [int(x) for x in user_args]
        ncsm_single_calculation(z=z0, a=a0, aeff=aeff0,
                                nhw=nhw0)
    elif len(user_args) == 5:
        z0, a0, aeff0, n1_0, n2_0 = [int(x) for x in user_args]
        ncsm_single_calculation(z=z0, a=a0, aeff=aeff0,
                                n1=n1_0, n2=n2_0)
    elif len(user_args) == 6:
        z0, a0, aeff0, nhw0, n1_0, n2_0 = [int(x) for x in user_args]
        ncsm_single_calculation(z=z0, a=a0, aeff=aeff0,
                                nhw=nhw0, n1=n1_0, n2=n2_0)
    else:
        raise InvalidNumberOfArgumentsException(
            '%d' % len(user_args) +
            ' is not a valid number of arguments for ncsm_single_calc.py.' +
            ' Please enter 2-6 arguments.'
        )
