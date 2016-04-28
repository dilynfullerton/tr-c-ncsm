import re
from os import path, listdir, remove, link, symlink, getcwd, makedirs
from FGetSmallerInteraction import run as truncate_interaction
from F_scale_interaction import scale_off_diag_outside_valence as scale_int


# CONSTANTS
MAX_NMAX = 15

# Directories
DPATH_MAIN = getcwd()
DPATH_TEMPLATES = path.join(DPATH_MAIN, 'templates')
DPATH_RESULTS = path.join(DPATH_MAIN, 'results')
DNAME_FMT_NUC = '%s%d_%d_Nhw%d_%d_%d'  # name, A, Aeff, nhw, n1, n2
DNAME_FMT_NUC_SF = DNAME_FMT_NUC + '_scale%.2f'  # scale factor
DNAME_FMT_VCE = 'vce_presc%d,%d,%d_Nmax%d_%d_%d_shell%d_dim%d'
#     A presc, Nmax, n1, n2, nshell, ncomponent
DNAME_FMT_VCE_SF = DNAME_FMT_VCE + '_scale%.2f'  # scale factor
DNAME_NCSD = 'ncsd'
DNAME_VCE = 'vce'

# Files
RGX_FNAME_TBME = 'TBME'
RGX_FNAME_TMP = '.*\.tmp'
FNAME_FMT_TBME = 'TBMEA2srg-n3lo2.O_%d.24'  # n1
FNAME_FMT_TBME_SF = FNAME_FMT_TBME + '_sf%.2f'  # n1 scalefactor
FNAME_FMT_NCSD_OUT = '%s%d_%d_Nhw%d_%d_%d.out'  # name, A, Aeff, Nhw, n1, n2
FNAME_FMT_NCSD_OUT_SF = FNAME_FMT_NCSD_OUT[:-4] + '_scale%.2f' + '.out'
FNAME_FMT_JOBSUB = FNAME_FMT_NCSD_OUT[:-4] + '.sh'
FNAME_FMT_JOBSUB_SF = FNAME_FMT_NCSD_OUT_SF[:-4] + '.sh'
FNAME_FMT_EGV = 'mfdp_%d.egv'  # Nhw
FNAME_TMP_MFDP = 'mfdp.dat'
FNAME_TMP_TRDENS_IN = 'trdens.in'
FNAME_TMP_JOBSUB = 'job.sh'
FNAME_TRDENS_IN = 'trdens.in'
LNAME_TBME = 'TBME.int'
LNAME_EGV = 'mfdp.egv'
LINE_FMT_MFDP_RESTR = ' %d %-2d %d %-2d %d %-4d ! N=%d'

# other
Z_NAME_MAP = {
    1: 'h_', 2: 'he', 3: 'li', 4: 'be', 5: 'b_', 6: 'c_', 7: 'n_', 8: 'o_',
    9: 'f_', 10: 'ne', 11: 'na', 12: 'mg', 13: 'al', 14: 'si', 15: 'p_',
    16: 's_', 17: 'cl', 18: 'ar', 19: 'k_', 20: 'ca'
}
Z_NAME_FMT_ALT = '%d-'
NCSD_NUM_STATES = 20
NCSD_NUM_ITER = 200


# FUNCTIONS
def _get_name(z, z_name_map=Z_NAME_MAP, alt_name=Z_NAME_FMT_ALT):
    """Given the proton number, return the short element name
    :param z: number of protons
    :param z_name_map: map from proton number to abbreviated name
    :param alt_name: alternate name format if z not in z_name_map
    """
    if z in z_name_map:
        return z_name_map[z]
    else:
        return alt_name % z


def _make_base_directories(a_aeff_to_dpath_map):
    """Makes directories for (A, Aeff) pairs if they do not exist yet
    :param a_aeff_to_dpath_map: map from A value to directory path
    directories are put
    """
    for d in a_aeff_to_dpath_map.values():
        if not path.exists(d):
            makedirs(d)


def make_vce_directories(
        a_prescription, nmax, n1, n2, nshell, ncomponent,
        dpath_results=DPATH_RESULTS, dname_vce=DNAME_VCE, scalefactor=None
):
    fmt = tuple(tuple(a_prescription) + (nmax, n1, n2, nshell, ncomponent))
    if scalefactor is not None:
        fmt += (scalefactor,)
        dname_fmt_vce = DNAME_FMT_VCE_SF
    else:
        dname_fmt_vce = DNAME_FMT_VCE
    vce_dirpath = path.join(dpath_results, dname_vce, dname_fmt_vce % fmt)
    if not path.exists(vce_dirpath):
        makedirs(vce_dirpath)
    return vce_dirpath


def _get_mfdp_restrictions_lines(nmax, max_allowed_nmax=MAX_NMAX):
    lines = list()
    for n in range(min(nmax, max_allowed_nmax) + 1):
        i = (n + 1) * (n + 2)
        lines.append(LINE_FMT_MFDP_RESTR % (0, i, 0, i, 0, 2 * i, n))
    return '\n'.join(lines)


class UnknownParityException(Exception):
    pass


def _get_parity(a, nshell):
    if nshell == 1:
        return a % 2
    elif nshell == 2:
        return 0
    else:
        raise UnknownParityException(
            '\nMethod of determining parity is not known for '
            'Nshell = %d' % nshell
        )


def _get_mfdp_replace_map(
        fname_tbme, outfile_name, z, a, n_hw, n_1, n_2, aeff, nshell,
        num_states=NCSD_NUM_STATES, num_iter=NCSD_NUM_ITER
):
    n = a - z
    par = _get_parity(a, nshell)
    tot2 = a % 2
    rest_lines = _get_mfdp_restrictions_lines(nmax=max(n_1, n_2))
    return {
        '<<TBMEFILE>>': str(fname_tbme),
        '<<OUTFILE>>': str(outfile_name),
        '<<Z>>': str(z), '<<N>>': str(n),
        '<<NHW>>': str(n_hw), '<<PAR>>': str(par), '<<TOT2>>': str(tot2),
        '<<N1>>': str(n_1), '<<N2>>': str(n_2),
        '<<RESTRICTIONS>>': str(rest_lines),
        '<<NUMST>>': str(num_states),
        '<<NUMITER>>': str(num_iter),
        '<<AEFF>>': str(aeff)
    }


def _rewrite_file(src, dst, replace_map):
    """Reads the file given by src, replaces string elements based
       on the replace map, writes the file into dst.
    """
    # read the src file
    fin = open(src, 'r')
    read_lines = fin.readlines()
    fin.close()
    # replace strings
    write_lines = list()
    for line in read_lines:
        for k, v in replace_map.iteritems():
            if k in line:
                line = line.replace(k, str(v))
        write_lines.append(line)
    # write to the dst file
    fout = open(dst, 'w')
    fout.writelines(write_lines)
    fout.close()


def _make_mfdp_file(
        z, a, aeff, nhw, nshell, n1, n2, dpath_elt, fname_outfile, fname_tbme,
        dpath_temp=DPATH_TEMPLATES,
        fname_tmp_mfdp=FNAME_TMP_MFDP,
):
    """Reads the mfdp file from path_temp
    and rewrites it into path_elt in accordance
    ith the given z, a, aeff, nhw, n1, n2, and outfile name
    :param z: Proton number
    :param a: Mass number
    :param aeff: Effective mass number for interaction
    :param nhw: Something something something dark side
    :param nshell: major oscillator shell (0=s, 1=p, 2=sd, ...)
    :param n1: Number of allowed states for single particles
    :param n2: Number of allowed states for two particles
    :param dpath_elt: path to the directory into which the mfdp file is
    to be put
    :param fname_outfile: name of the outfile
    :param dpath_temp: path to the template directory
    :param fname_tmp_mfdp: name of the mfdp file
    """
    temp_mfdp_path = path.join(dpath_temp, fname_tmp_mfdp)
    mfdp_path = path.join(dpath_elt, fname_tmp_mfdp)
    replace_map = _get_mfdp_replace_map(
        fname_tbme=fname_tbme, outfile_name=fname_outfile,
        z=z, a=a, n_hw=nhw, nshell=nshell, n_1=n1, n_2=n2, aeff=aeff
    )
    _rewrite_file(src=temp_mfdp_path, dst=mfdp_path, replace_map=replace_map)


def _make_mfdp_files(
        a_list, aeff_list, nhw_list, z, nshell, n_1, n_2, fname_tbme,
        a_aeff_to_dpath_map, a_aeff_to_outfile_fpath_map,
        dpath_temp=DPATH_TEMPLATES, fname_tmp_mfdp=FNAME_TMP_MFDP,
):
    for a, aeff, nhw in zip(a_list, aeff_list, nhw_list):
        outfile_name = path.split(a_aeff_to_outfile_fpath_map[(a, aeff)])[1]
        _make_mfdp_file(
            z=z, a=a, aeff=aeff, nhw=nhw, nshell=nshell, n1=n_1, n2=n_2,
            dpath_elt=a_aeff_to_dpath_map[(a, aeff)], dpath_temp=dpath_temp,
            fname_outfile=outfile_name, fname_tmp_mfdp=fname_tmp_mfdp,
            fname_tbme=fname_tbme,
        )


class UnknownNumStatesException(Exception):
    pass


# todo IMPORTANT write this method to be general
# This should determine 'nhw_mod' and 'dim_nhw_mod'. These are the names
# of the variables in trdens-kernels.f, which are retrieved from trdens.in.
# I do not know what they mean, and the means of determining them
# is a hard-coded case-by-case system that throws an exception if it does
# not know what to do. Also, it could potentially be incorrect and not
# throw an exception if more variables are needed to determine the result.
# This REALLY needs fixing.
def _get_num_states(z, a, a0):
    nhw_mod = a - a0
    if nhw_mod == 0:
        dim_nhw_mod = 1
    elif nhw_mod == 1:
        dim_nhw_mod = 2
    elif nhw_mod == 2:
        dim_nhw_mod = 5
    else:
        raise UnknownNumStatesException()
    return nhw_mod, dim_nhw_mod


def _get_trdens_replace_map(z, a, a0):
    nnn, num_states = _get_num_states(z, a, a0)
    return {'<<NNN>>': str(nnn), '<<NUMSTATES>>': str(num_states)}


def make_trdens_file(
        z, a, a0, nuc_dir,
        dpath_results=DPATH_RESULTS, dpath_temp=DPATH_TEMPLATES,
        fname_tmp_trdens_in=FNAME_TMP_TRDENS_IN
):
    """Reads the trdens.in file from path_temp and rewrites it
    into path_elt in accordance with the given z, a
    :param z: proton number
    :param a: mass number
    :param a0: lowest mass number in the shell
    :param nuc_dir: directory name
    :param dpath_results: path to the results directory
    :param dpath_temp: path to the templates directory
    :param fname_tmp_trdens_in: name of the trdens file in the templates dir
    """
    src = path.join(dpath_temp, fname_tmp_trdens_in)
    path_elt = path.join(dpath_results, nuc_dir)
    dst = path.join(path_elt, FNAME_TRDENS_IN)
    rep_map = _get_trdens_replace_map(z=z, a=a, a0=a0)
    _rewrite_file(src=src, dst=dst, replace_map=rep_map)


class TbmeFileNotFoundException(Exception):
    pass


def _truncate_space(
        nshell, n1, n2, dpath_elt, scalefactor,
        dpath_templates=DPATH_TEMPLATES, force=False,
):
    """Run the script that truncates the space by removing extraneous
    interactions from the TBME file
    :param n1: Maximum state for single particle
    :param n2: Maximum state for two particles
    :param dpath_elt: Path to the directory in which the resultant TBME file
    is to be put
    """
    filenames = listdir(dpath_templates)
    for f in filenames:
        if re.match(RGX_FNAME_TBME, f):
            tbme_filename = f
            break
    else:
        raise TbmeFileNotFoundException(
            '\nTBME interaction file not found in directory '
            '%s' % dpath_templates
        )
    src_path = path.join(dpath_templates, tbme_filename)
    tmp_path = path.join(dpath_elt, FNAME_FMT_TBME % n1)
    link_path = path.join(dpath_elt, LNAME_TBME)
    if path.exists(link_path):
        remove(link_path)
    if scalefactor is None:
        dst_path = tmp_path
    else:
        dst_path = path.join(dpath_elt, FNAME_FMT_TBME_SF % (n1, scalefactor))
    if force or not path.exists(dst_path):
        if path.exists(dst_path):
            remove(dst_path)
        truncate_interaction(src_path, n1, n2, tmp_path)
        if scalefactor is not None:
            scale_int(src=tmp_path, dst=dst_path, nshell=nshell,
                      scalefactor=scalefactor)
            remove(tmp_path)
    symlink(dst_path, link_path)
    return dst_path, link_path


def _get_job_replace_map(walltime):
    return {'<<WALLTIME>>': str(walltime)}


def _make_job_submit_file(
        dst_fpath, walltime,
        dpath_temp=DPATH_TEMPLATES, fname_tmp_jobsub=FNAME_TMP_JOBSUB
):
    src_fpath = path.join(dpath_temp, fname_tmp_jobsub)
    rep_map = _get_job_replace_map(walltime=walltime)
    _rewrite_file(src=src_fpath, dst=dst_fpath, replace_map=rep_map)


def _make_job_submit_files(a_aeff_to_jobsub_fpath_map, walltime):
    job_files = list(a_aeff_to_jobsub_fpath_map.values())
    dst = job_files.pop()
    _make_job_submit_file(dst_fpath=dst, walltime=walltime)
    # link this file to rest of destination paths
    for dst2 in job_files:
        if path.exists(dst2):
            remove(dst2)
        link(dst, dst2)


def _truncate_spaces(nshell, n1, n2, dirpaths, scalefactor, force=False):
    """For multiple directories, perform the operation of truncate_space
    :param n1: max allowed one-particle state
    :param n2: max allowed two-particle state
    :param dirpaths: Paths to the destination directories
    """
    dirpaths = list(dirpaths)
    d0 = dirpaths.pop()
    # truncate interaction once
    fpath0, lpath0 = _truncate_space(
        nshell=nshell, n1=n1, n2=n2, dpath_elt=d0, scalefactor=scalefactor,
        force=force
    )
    fname_tbme = path.split(fpath0)[1]
    lname_tbme = path.split(lpath0)[1]
    # link truncated interaction file to the rest of the directories
    for d in dirpaths:
        dst_path = path.join(d, fname_tbme)
        sl_path = path.join(d, lname_tbme)
        if path.exists(sl_path):  # symlink exists
            remove(sl_path)
        if not path.exists(dst_path):  # link to TBME does not exist
            link(fpath0, dst_path)
        elif force:
            remove(dst_path)
            link(fpath0, dst_path)
        symlink(dst_path, sl_path)
    return fname_tbme, lname_tbme


class EgvFileNotFoundException(Exception):
    pass


def rename_egv_file(a6_dir, nhw, force):
    """Renames the egv file from its default output name to the name needed
    for running TRDENS
    :param a6_dir: Directory in which the file resides
    :param nhw: major oscillator model space truncation
    :param force: If True, replaces any existing files by the name of
    next_egv_name
    """
    next_egv_path = path.join(a6_dir, LNAME_EGV)
    if path.lexists(next_egv_path):
        if not force:
            return 0
        else:
            remove(next_egv_path)
    filenames = listdir(a6_dir)
    fname_egv = FNAME_FMT_EGV % nhw
    for f in filenames:
        if f == fname_egv:
            symlink(f, next_egv_path)
            break
    else:
        raise EgvFileNotFoundException('File not found: %s' % fname_egv)
    return 1


def get_a_aeff_to_dpath_map(
        a_list, aeff_list, nhw_list, z, n1, n2,
        dpath_results=DPATH_RESULTS, dname_ncsd=DNAME_NCSD, scalefactor=None,
):
    a_paths_map = dict()
    if scalefactor is not None:
        dname_fmt_nuc = DNAME_FMT_NUC_SF
    else:
        dname_fmt_nuc = DNAME_FMT_NUC
    path_fmt = path.join(dpath_results, dname_ncsd, dname_fmt_nuc)
    for a, aeff, nhw in zip(a_list, aeff_list, nhw_list):
        args = (_get_name(z), a, aeff, nhw, n1, n2)
        if scalefactor is not None:
            args += (scalefactor,)
        a_paths_map[(a, aeff)] = path_fmt % args
    return a_paths_map


def get_a_aeff_to_outfile_fpath_map(
        a_list, aeff_list, nhw_list, z, n1, n2, a_aeff_to_dirpath_map,
        scalefactor=None
):
    a_aeff_outfile_map = dict()
    if scalefactor is not None:
        fname_fmt = FNAME_FMT_NCSD_OUT_SF
    else:
        fname_fmt = FNAME_FMT_NCSD_OUT
    for a, aeff, nhw in zip(a_list, aeff_list, nhw_list):
        args = (_get_name(z), a, aeff, nhw, n1, n2)
        if scalefactor is not None:
            args += (scalefactor,)
        a_aeff_outfile_map[(a, aeff)] = path.join(
            a_aeff_to_dirpath_map[(a, aeff)], fname_fmt % args)
    return a_aeff_outfile_map


def _get_a_aeff_to_egv_fpath_map(
        a_list, aeff_list, nhw_list, a_aeff_to_dirpath_map):
    a_aeff_to_egv_map = dict()
    for a, aeff, nhw in zip(a_list, aeff_list, nhw_list):
        a_aeff_to_egv_map[(a, aeff)] = path.join(
            a_aeff_to_dirpath_map[(a, aeff)], FNAME_FMT_EGV % nhw)
    return a_aeff_to_egv_map


def _get_a_aeff_to_jobsub_fpath_map(
        a_list, aeff_list, nhw_list, z, n1, n2, a_aeff_to_dirpath_map,
        scalefactor=None,
):
    a_aeff_jobsub_map = dict()
    if scalefactor is not None:
        fname_fmt = FNAME_FMT_JOBSUB_SF
    else:
        fname_fmt = FNAME_FMT_JOBSUB
    for a, aeff, nhw in zip(a_list, aeff_list, nhw_list):
        args = (_get_name(z), a, aeff, nhw, n1, n2)
        if scalefactor is not None:
            args += (scalefactor,)
        a_aeff_jobsub_map[(a, aeff)] = path.join(
            a_aeff_to_dirpath_map[(a, aeff)], fname_fmt % args)
    return a_aeff_jobsub_map


def remove_ncsd_tmp_files(dpaths_list):
    """Removes all *.tmp files from the given list of directories
    :param dpaths_list: list of absolute directory paths
    """
    for dpath in dpaths_list:
        for fname in listdir(dpath):
            if re.match(RGX_FNAME_TMP, fname):
                remove(path.join(dpath, fname))


def prepare_directories(
        a_list, aeff_list, nhw_list, z, n1, n2, nshell, scalefactor,
        dpath_templates, dpath_results,
        cluster_submit=False, walltime=None, progress=False,
        force=False,
):
    """Creates directories and files necessary to run NCSD calculations.
    Returns maps to the important files and directories.
    :param a_list: ordered list of A values for which to create directories
    :param aeff_list: ordered list of Aeff values (corresponding to the
    respective A values in a_list) for which to create directories
    :param nhw_list: ordered list of Nhw values (corresponding to the
    respective A and Aeff values in the a_list and aeff_list) for which to
    create directories
    :param z: proton number Z
    :param n1: max allowed 1-particle state
    :param n2: max allowed 2-particle state
    :param nshell: shell (0: s, 1: p, 2: sd, ...)
    :param scalefactor: factor by which to scale off-diagonal coupling terms of
    the TBME interaction
    :param dpath_templates: path to the templates directory
    :param dpath_results: path to the results directory
    :param progress: if true, shows a progress bar (verbose mode will be off)
    :param cluster_submit: if true, submits the job to the OpenMP cluster
    using qsub
    :param walltime: walltime to be allotted to a cluster submission
    :param force: if true, forces re-truncation of the TBME interaction file
    :return (A,Aeff)->dir, (A,Aeff)->*.egv, (A,Aeff)->*.sh, (A,Aeff)->*.out
    """
    # get maps
    a_aeff_to_dir_map = get_a_aeff_to_dpath_map(
        a_list=a_list, aeff_list=aeff_list, nhw_list=nhw_list,
        z=z, n1=n1, n2=n2, scalefactor=scalefactor,
        dpath_results=dpath_results,
    )
    a_aeff_to_outfile_map = get_a_aeff_to_outfile_fpath_map(
        a_list=a_list, aeff_list=aeff_list, nhw_list=nhw_list,
        z=z, n1=n1, n2=n2, a_aeff_to_dirpath_map=a_aeff_to_dir_map,
        scalefactor=scalefactor,
    )
    a_aeff_to_egvfile_map = _get_a_aeff_to_egv_fpath_map(
        a_list=a_list, aeff_list=aeff_list, nhw_list=nhw_list,
        a_aeff_to_dirpath_map=a_aeff_to_dir_map
    )
    if cluster_submit:
        a_aeff_to_jobfile_map = _get_a_aeff_to_jobsub_fpath_map(
            a_list=a_list, aeff_list=aeff_list, nhw_list=nhw_list,
            z=z, n1=n1, n2=n2, a_aeff_to_dirpath_map=a_aeff_to_dir_map,
            scalefactor=scalefactor,
        )
    else:
        a_aeff_to_jobfile_map = dict()
    # make stuff
    if progress:
        print '  Making directories...'
    _make_base_directories(a_aeff_to_dpath_map=a_aeff_to_dir_map)
    if progress:
        print '  Truncating interaction to N1=%d N2=%d...' % (n1, n2)
    fname_tbme, lname_tbme = _truncate_spaces(
        nshell=nshell, n1=n1, n2=n2, dirpaths=a_aeff_to_dir_map.values(),
        scalefactor=scalefactor, force=force
    )
    if progress:
        print '  Writing mfdp files...'
    _make_mfdp_files(
        a_list=a_list, aeff_list=aeff_list, nhw_list=nhw_list,
        a_aeff_to_dpath_map=a_aeff_to_dir_map,
        a_aeff_to_outfile_fpath_map=a_aeff_to_outfile_map,
        fname_tbme=lname_tbme, dpath_temp=dpath_templates,
        z=z, nshell=nshell, n_1=n1, n_2=n2,
    )
    if cluster_submit:
        if progress:
            print '  Writing cluster submit files...'
        _make_job_submit_files(
            a_aeff_to_jobsub_fpath_map=a_aeff_to_jobfile_map,
            walltime=walltime,
        )
    return (a_aeff_to_dir_map, a_aeff_to_egvfile_map, a_aeff_to_jobfile_map,
            a_aeff_to_outfile_map)