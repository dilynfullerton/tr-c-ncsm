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
DNAME_FMT_NUC = '%s%d_%d_Nhw%d_%d_%d'
# name, A, Aeff, nhw, n1, n2,
DNAME_ADD_NUC_IPROT = '_ip0'
DNAME_ADD_NUC_SF = '_scale%.2f'  # scale factor
DNAME_FMT_VCE = 'vce_presc%d,%d,%d_Nmax%d_%d_%d_shell%d_dim%d'
#     A presc, Nmax, n1, n2, nshell, ncomponent
DNAME_ADD_VCE_IPROT = '_ip0'
DNAME_ADD_VCE_SF = '_scale%.2f'  # scale factor
DNAME_NCSD = 'ncsd'
DNAME_VCE = 'vce'

# Files
RGX_FNAME_TBME = 'TBME'
RGX_FNAME_TMP = '.*\.tmp'
FNAME_FMT_TBME = 'TBMEA2srg-n3lo2.O_%d.24'  # n1
FNAME_ADD_TBME_IPROT = '_ip0'
FNAME_ADD_TBME_SF = '_sf%.2f'  # scalefactor
FNAME_FMT_NCSD_OUT = DNAME_FMT_NUC
# name, A, Aeff, Nhw, n1, n2,
FNAME_ADD_NCSD_OUT_IPROT = DNAME_ADD_NUC_IPROT
FNAME_ADD_NCSD_OUT_SF = DNAME_ADD_NUC_SF
FNAME_EXT_NCSD_OUT = '.out'
FNAME_FMT_JOBSUB = FNAME_FMT_NCSD_OUT
FNAME_ADD_JOBSUB_IPROT = FNAME_ADD_NCSD_OUT_IPROT
FNAME_ADD_JOBSUB_SF = FNAME_ADD_NCSD_OUT_SF
FNAME_EXT_JOBSUB = '.sh'
FNAME_FMT_EGV = 'mfdp_%d.egv'  # Nhw
FNAME_MFDP = 'mfdp.dat'
FNAME_TMP_TRDENS_IN = 'trdens.in'
FNAME_TMP_JOBSUB_NCSD = 'ncsd_job.sh'
FNAME_TMP_JOBSUB_TRDENS = 'trdens_job.sh'
FNAME_TRDENS_IN = 'trdens.in'
LNAME_TBME = 'TBME.int'
LNAME_EGV = 'mfdp.egv'
LINE_FMT_MFDP_RESTR = ' %d %-2d %d %-2d %d %-4d ! N=%d'

# other
Z_NAME_MAP = {
    1: 'h-', 2: 'he', 3: 'li', 4: 'be', 5: 'b-', 6: 'c-', 7: 'n-', 8: 'o-',
    9: 'f-', 10: 'ne', 11: 'na', 12: 'mg', 13: 'al', 14: 'si', 15: 'p-',
    16: 's-', 17: 'cl', 18: 'ar', 19: 'k-', 20: 'ca'
}
Z_NAME_FMT_ALT = '%d-'


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


# TODO: does this really need its own function?
def _make_base_directories(a_aeff_to_dpath_map):
    """Makes directories for (A, Aeff) pairs if they do not exist yet
    :param a_aeff_to_dpath_map: map from A value to directory path
    directories are put
    """
    for d in a_aeff_to_dpath_map.values():
        if not path.exists(d):
            makedirs(d)


def make_vce_directory(
        a_prescription, nmax, n1, n2, nshell, ncomponent,
        dpath_results=DPATH_RESULTS, dname_vce=DNAME_VCE, scalefactor=None,
        remove_protons=False,
):
    """Makes the VCE directory for the given A-prescription
    :param a_prescription: 3-tuples of effective mass numbers, which are used
    in place of the first 3 mass numbers in the shell to generate the
    effective interaction
    :param nmax: oscillator model space truncation
    :param nshell: shell number (0 = s, 1 = p, 2 = sd, ...)
    :param ncomponent: dimension (1 = neutrons, 2 = protons and neutrons)
    :param n1: max allowed 1-particle state (for TBME truncation)
    :param n2: max allowed 2-particle state (for TBME truncation)
    :param dpath_results: path to the results directory
    :param dname_vce: name of the vce subdirectory
    :param scalefactor: factor by which off-diagonal coupling terms in the
    TBME interaction were scaled
    :param remove_protons: boolean indicating whether or not the proton part
    of the TBME interaction file was scaled to 0
    :return: path to the particular VCE directory in which interaction
    files for this prescription are to be saved
    """
    fmt = tuple(tuple(a_prescription) + (nmax, n1, n2, nshell, ncomponent))
    dname_fmt_vce = DNAME_FMT_VCE
    if remove_protons:
        dname_fmt_vce += DNAME_ADD_VCE_IPROT
    if scalefactor is not None:
        dname_fmt_vce += DNAME_ADD_VCE_SF
        fmt += (scalefactor,)
    vce_dirpath = path.join(dpath_results, dname_vce, dname_fmt_vce % fmt)
    if not path.exists(vce_dirpath):
        makedirs(vce_dirpath)
    return vce_dirpath


def _get_mfdp_restrictions_lines(nmax, max_allowed_nmax=MAX_NMAX):
    """Returns the mfdp lines, which specify occupation restriction.
    :param nmax: major oscillator model space truncation
    :param max_allowed_nmax: max number of ines to be generated
    """
    # TODO: Is this correct to produce restrictions up to Nmax, or should they
    # TODO: be produces up to Nhw?
    lines = list()
    for n in range(min(nmax, max_allowed_nmax) + 1):
        i = (n + 1) * (n + 2)
        lines.append(LINE_FMT_MFDP_RESTR % (0, i, 0, i, 0, 2 * i, n))
    return '\n'.join(lines)


class UnknownParityException(Exception):
    pass


# TODO: make this function general
def _get_parity(a, nshell):
    """Returns the parity for a particular mass number and shell
    :param a: mass number (A)
    :param nshell: major oscillator shell (0=s, 1=p, 2=sd, ...)
    """
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
        fname_tbme, outfile_name, z, a, n_hw, n_1, n_2, aeff, nshell, beta_cm,
        num_states, num_iter,
):
    """Returns a map from the placeholder string in mfdp.dat template file to
    the value or string that should replace it.
    :param fname_tbme: name of the TBME interaction file
    :param outfile_name: name of the NCSD *.out file
    :param z: proton number (Z)
    :param a: mass number (A)
    :param n_hw: major oscillator truncation Nhw
    :param n_1: max allowed 1-particle state (for TBME truncation)
    :param n_2: max allowed 2-particle state (for TBME truncation)
    :param aeff: effective mass number
    :param nshell: major oscillator shell (0=s, 1=p, 2=sd, ...)
    :param beta_cm: center of mass beta term
    :param num_states: max number of states to calculate
    :param num_iter: max number of iteractions for lanczos
    """
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
        '<<BETA>>': str(float(beta_cm)),
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
        z, a, aeff, nhw, nshell, n1, n2, beta_cm, num_states, num_iter,
        dpath_elt, fname_outfile, fname_tbme, dpath_temp,
):
    """Reads the mfdp file from path_temp
    and rewrites it into path_elt in accordance
    ith the given z, a, aeff, nhw, n1, n2, and outfile name
    :param z: Proton number (Z)
    :param a: Mass number (A)
    :param aeff: Effective mass number for interaction
    :param nhw: major oscillator truncation
    :param nshell: major oscillator shell (0=s, 1=p, 2=sd, ...)
    :param n1: Number of allowed states for single particles
    :param n2: Number of allowed states for two particles
    :param beta_cm: center of mass beta parameter
    :param num_states: number of NCSD states to calculate
    :param num_iter: number of iteractions for lanczos algorithm
    :param dpath_elt: path to the directory into which the mfdp file is
    to be put
    :param fname_outfile: name of the outfile
    :param dpath_temp: path to the template directory
    """
    temp_mfdp_path = path.join(dpath_temp, FNAME_MFDP)
    mfdp_path = path.join(dpath_elt, FNAME_MFDP)
    replace_map = _get_mfdp_replace_map(
        fname_tbme=fname_tbme, outfile_name=fname_outfile, z=z, a=a, n_hw=nhw,
        nshell=nshell, n_1=n1, n_2=n2, aeff=aeff, beta_cm=beta_cm,
        num_states=num_states, num_iter=num_iter,
    )
    _rewrite_file(src=temp_mfdp_path, dst=mfdp_path, replace_map=replace_map)


def _make_mfdp_files(
        a_list, aeff_list, nhw_list, z, nshell, n_1, n_2, beta_cm, num_states,
        num_iter, a_aeff_to_dpath_map, a_aeff_to_outfile_fpath_map, fname_tbme,
        dpath_temp,
):
    """Writes mfdp files for the specified NCSD calculations
    :param a_list: ordered list of A values
    :param aeff_list: ordered list of Aeff values, which match to the a_list
    :param nhw_list: ordered list of Nhw values, which match to the a_list
    :param z: proton number (Z)
    :param nshell: major oscillator shell (0=s, 1=p, 2=sd,...)
    :param n_1: max allowed 1-particle state (for TBME truncation)
    :param n_2: max allowed 2-particle state (for TBME truncation)
    :param beta_cm: ceter of mass beta term
    :param num_states: max number of NCSD states to calculate
    :param num_iter: number of iterations for lanczos algorithm
    :param a_aeff_to_dpath_map: map from (A,Aeff) to the directory in which
    the NCSD calculation was done
    :param a_aeff_to_outfile_fpath_map: map from (A,Aeff) to the *.out file
    to be generated by NCSD
    :param fname_tbme: name of the TBME interaction file
    :param dpath_temp: path to the templates directory in which the mfdp.dat
    template file is stored
    """
    for a, aeff, nhw in zip(a_list, aeff_list, nhw_list):
        outfile_name = path.split(a_aeff_to_outfile_fpath_map[(a, aeff)])[1]
        _make_mfdp_file(
            z=z, a=a, aeff=aeff, nhw=nhw, nshell=nshell, n1=n_1, n2=n_2,
            beta_cm=beta_cm, num_states=num_states, num_iter=num_iter,
            dpath_elt=a_aeff_to_dpath_map[(a, aeff)], dpath_temp=dpath_temp,
            fname_outfile=outfile_name, fname_tbme=fname_tbme,
        )


class UnknownNumStatesException(Exception):
    pass


# TODO: IMPORTANT write this method to be general
# This should determine 'nhw_mod' and 'dim_nhw_mod'. These are the names
# of the variables in trdens-kernels.f, which are retrieved from trdens.in.
# I am only somewhat confident that this works when A - A0 = 0
# (he6 in p shell, o18 in sd shell, etc). For anything else, I would not trust
# it.
def _get_num_states(a, a0, nshell, nmax):
    """Given a mass number, the first mass number in the shell, and the shell,
    returns the model space Nhw and model space dimension for use in the
    trdens.in input file
    :param a: mass number (A)
    :param a0: first mass in the shell
    :param nshell: major oscillator shell (0=s, 1=p, 2=sd,...)
    :param nmax: shell truncation
    """
    nhw_mod = (a - a0) * nshell  # + nmax TODO: is this correct generally
    dim_nhw_mod = 0
    for j1_2 in range(1, nhw_mod+2, 2):
        dim_nhw_mod += (j1_2 + 1)/2
        for j2_2 in range(j1_2+1, nhw_mod+2, 2):
            dim_nhw_mod += j1_2 + 1
    return nhw_mod, dim_nhw_mod


def _get_trdens_replace_map(a, a0, nshell, nmax):
    """Given a mass number, the first mass number in the shell, and the shell,
    returns a map from placeholder string in the trdens.in template file to
    the value or string that replaces it.
    :param a: mass number (A)
    :param a0: first mass in the shell
    :param nshell: major oscillator shell (0=s, 1=p, 2=sd,...)
    :param nmax: shell truncation
    """
    nnn, num_states = _get_num_states(a, a0, nshell, nmax)
    return {'<<NNN>>': str(nnn), '<<NUMSTATES>>': str(num_states)}


def make_trdens_file(
        a, a0, nuc_dir, nshell, nmax, dpath_results=DPATH_RESULTS,
        dpath_temp=DPATH_TEMPLATES, fname_tmp_trdens_in=FNAME_TMP_TRDENS_IN
):
    """Reads the trdens.in file from path_temp and rewrites it
    into path_elt in accordance with the given z, a
    :param a: mass number
    :param a0: lowest mass number in the shell
    :param nshell: major oscillator shell (0=s, 1=p, 2=sd, ...)
    :param nmax: shell truncation
    :param nuc_dir: directory name
    :param dpath_results: path to the results directory
    :param dpath_temp: path to the templates directory
    :param fname_tmp_trdens_in: name of the trdens file in the templates dir
    """
    src = path.join(dpath_temp, fname_tmp_trdens_in)
    path_elt = path.join(dpath_results, nuc_dir)
    dst = path.join(path_elt, FNAME_TRDENS_IN)
    rep_map = _get_trdens_replace_map(a=a, a0=a0, nshell=nshell, nmax=nmax)
    _rewrite_file(src=src, dst=dst, replace_map=rep_map)


class TbmeFileNotFoundException(Exception):
    pass


def _truncate_space(
        nshell, n1, n2, dpath_elt, scalefactor, remove_protons,
        dpath_templates=DPATH_TEMPLATES, force=False,
):
    """Run the script that truncates the space by removing extraneous
    interactions from the TBME file
    :param nshell: major oscillator shell (0=s, 1=p, 2=sd, 3=fp, ...)
    :param n1: Maximum state for single particle
    :param n2: Maximum state for two particles
    :param dpath_elt: Path to the directory in which the resultant TBME file
    is to be put
    :raises TbmeFileNotFoundException: when the TBME file is not found in
    the templates directory
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
    # Determine the final destination path
    dst_path = tmp_path
    if remove_protons:
        dst_path += FNAME_ADD_TBME_IPROT
    if scalefactor is not None:
        dst_path += FNAME_ADD_TBME_SF % scalefactor
    # Remove symlink if already exists
    if path.lexists(link_path) or path.exists(link_path):
        remove(link_path)
    if path.split(link_path)[1] in listdir(dpath_elt):
        remove(link_path)
    # Do truncation [and scaling]
    if force or not path.exists(dst_path):
        if path.exists(dst_path):
            remove(dst_path)
        if n1 != 15 or n2 != 15:
            truncate_interaction(src_path, n1, n2, tmp_path)
        else:
            symlink(src_path, tmp_path)
        if (scalefactor is not None and scalefactor != 1) or remove_protons:
            if scalefactor is None:
                scalefactor = 1.0
            scale_int(
                src=tmp_path, dst=dst_path, nshell=nshell,
                scalefactor=scalefactor, rm_proton_interaction=remove_protons
            )
            remove(tmp_path)
    symlink(dst_path, link_path)
    return dst_path, link_path


def _get_job_replace_map(walltime):
    """Returns a map from placeholder string in the job.sh template file to
    the string that is to replace it
    :param walltime: string in format hh:mm:ss that represents the amount of
    walltime a job is to be given when submitted to the cluster
    """
    return {'<<WALLTIME>>': str(walltime)}


def _make_job_submit_file(
        dst_fpath, walltime, fname_tmp_jobsub, dpath_temp=DPATH_TEMPLATES,
):
    """For the given walltime, writes the job.sh submit file to the given
    destination (dst_fpath)
    :param dst_fpath: destination file path
    :param walltime: string representation of the allotted time on the cluster
    :param dpath_temp: directory in which the job.sh template file is stored
    :param fname_tmp_jobsub: name of the job.sh template file
    """
    src_fpath = path.join(dpath_temp, fname_tmp_jobsub)
    rep_map = _get_job_replace_map(walltime=walltime)
    _rewrite_file(src=src_fpath, dst=dst_fpath, replace_map=rep_map)


def _make_job_submit_files(
        a_aeff_to_jobsub_fpath_map, walltime, fname_tmp_jobsub
):
    """For each jobsub_fpath in a_aeff_to_jobsub_fpath_map, makes the job.sh
    file for submission to the cluster. As these are all the same for a given
    walltime, one file is written and the rest are linked.
    :param a_aeff_to_jobsub_fpath_map: map from (A,Aeff) to the job submission
    file
    :param walltime: string representation of the allowed walltime in the
    format hh:mm:ss
    """
    job_files = list(a_aeff_to_jobsub_fpath_map.values())
    dst = job_files.pop()
    _make_job_submit_file(
        dst_fpath=dst, walltime=walltime, fname_tmp_jobsub=fname_tmp_jobsub)
    # link this file to rest of destination paths
    for dst2 in job_files:
        if path.exists(dst2):
            remove(dst2)
        link(dst, dst2)


def _truncate_spaces(
        nshell, n1, n2, dirpaths, scalefactor, remove_protons, force=False):
    """For multiple directories, perform the operation of truncate_space
    :param nshell: major oscillator shell (0=s, 1=p, ...)
    :param n1: max allowed one-particle state
    :param n2: max allowed two-particle state
    :param dirpaths: Paths to the destination directories
    """
    dirpaths = list(dirpaths)
    d0 = dirpaths.pop()
    # truncate interaction once
    fpath0, lpath0 = _truncate_space(
        nshell=nshell, n1=n1, n2=n2, dpath_elt=d0, scalefactor=scalefactor,
        remove_protons=remove_protons, force=force,
    )
    fname_tbme = path.split(fpath0)[1]
    lname_tbme = path.split(lpath0)[1]
    # link truncated interaction file to the rest of the directories
    for d in dirpaths:
        dst_path = path.join(d, fname_tbme)
        sl_path = path.join(d, lname_tbme)
        try:
            if path.exists(sl_path) or path.lexists(sl_path):  # symlink exists
                remove(sl_path)
            if not (path.exists(dst_path) or path.lexists(dst_path)):
                link(fpath0, dst_path)
            elif force:
                remove(dst_path)
                link(fpath0, dst_path)
        except OSError:
            print 'Could not link %s to %s.' % (fpath0, dst_path)
            raise
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
        a_list, aeff_list, nhw_list, z, n1, n2, scalefactor=None,
        remove_protons=False, dpath_results=DPATH_RESULTS,
        dname_ncsd=DNAME_NCSD,
):
    """Given ordered lists a_list, aeff_list, nhw_list, constructs a map from
    (A,Aeff) to the directory in which the NCSD calculation is to be done
    :param a_list: ordered list of mass numbers (A)
    :param aeff_list: ordered list of effective mass numbers (Aeff),
    corresponding to the mass numbers in a_list
    :param nhw_list: ordered list of Nhw corresponding to the mass numbers in
    a_list
    :param z: proton number (Z)
    :param n1: max allowed one-particle state (for TBME truncation)
    :param n2: max allowed two-particle state (for TBME truncation)
    :param dpath_results: path to the results directory
    :param dname_ncsd: name of the NCSD subdirectory
    :param scalefactor: value by which off-diagonal coupling terms in the
    TBME interaction were scaled.
    :param remove_protons: either true or false, indicating whether the
    proton parts of the interaction (Vpp and Vpn) are being scaled to 0
    """
    a_paths_map = dict()
    dname_fmt_nuc = DNAME_FMT_NUC
    if remove_protons:
        dname_fmt_nuc += DNAME_ADD_NUC_IPROT
    if scalefactor is not None:
        dname_fmt_nuc += DNAME_ADD_NUC_SF
    path_fmt = path.join(dpath_results, dname_ncsd, dname_fmt_nuc)
    for a, aeff, nhw in zip(a_list, aeff_list, nhw_list):
        args = (_get_name(z), a, aeff, nhw, n1, n2)
        if scalefactor is not None:
            args += (scalefactor,)
        a_paths_map[(a, aeff)] = path_fmt % args
    return a_paths_map


def get_a_aeff_to_outfile_fpath_map(
        a_list, aeff_list, nhw_list, z, n1, n2, remove_protons,
        a_aeff_to_dirpath_map, scalefactor=None
):
    """Given ordered lists a_list, aeff_list, nhw_list, constructs a map from
    (A,Aeff) to the NCSD *.out file that will be created
    :param a_list: ordered list of mass numbers (A)
    :param aeff_list: ordered list of effective mass numbers (Aeff),
    corresponding to the mass numbers in a_list
    :param nhw_list: ordered list of Nhw corresponding to the mass numbers in
    a_list
    :param z: proton number (Z)
    :param n1: max allowed one-particle state (for TBME truncation)
    :param n2: max allowed two-particle state (for TBME truncation)
    :param a_aeff_to_dirpath_map: map from (A,Aeff) to the directory in which
    an NCSD calculation will be done
    :param scalefactor: value by which off-diagonal coupling terms in the
    TBME interaction were scaled.
    :param remove_protons: either true or false, indicating whether the
    proton parts of the interaction (Vpp and Vpn) are being scaled to 0
    """
    a_aeff_outfile_map = dict()
    fname_fmt = FNAME_FMT_NCSD_OUT
    if remove_protons:
        fname_fmt += FNAME_ADD_NCSD_OUT_IPROT
    if scalefactor is not None:
        fname_fmt += FNAME_ADD_NCSD_OUT_SF
    fname_fmt += FNAME_EXT_NCSD_OUT
    for a, aeff, nhw in zip(a_list, aeff_list, nhw_list):
        args = (_get_name(z), a, aeff, nhw, n1, n2)
        if scalefactor is not None:
            args += (scalefactor,)
        a_aeff_outfile_map[(a, aeff)] = path.join(
            a_aeff_to_dirpath_map[(a, aeff)], fname_fmt % args)
    return a_aeff_outfile_map


def _get_a_aeff_to_egv_fpath_map(
        a_list, aeff_list, nhw_list, a_aeff_to_dirpath_map):
    """Given ordered lists a_list, aeff_list, nhw_list, constructs a map from
    (A,Aeff) to the *.egv file, created by NCSD, which is used to signify that
    the calculation has been completed
    :param a_list: ordered list of mass numbers (A)
    :param aeff_list: ordered list of effective mass numbers (Aeff),
    corresponding to the mass numbers in a_list
    :param nhw_list: ordered list of Nhw corresponding to the mass numbers in
    a_list
    :param a_aeff_to_dirpath_map: map from (A,Aeff) to the directory in which
    an NCSD calculation will be done
    """
    a_aeff_to_egv_map = dict()
    for a, aeff, nhw in zip(a_list, aeff_list, nhw_list):
        a_aeff_to_egv_map[(a, aeff)] = path.join(
            a_aeff_to_dirpath_map[(a, aeff)], FNAME_FMT_EGV % nhw)
    return a_aeff_to_egv_map


def _get_a_aeff_to_jobsub_fpath_map(
        a_list, aeff_list, nhw_list, z, n1, n2, remove_protons,
        a_aeff_to_dirpath_map, scalefactor=None,
):
    """Given ordered lists a_list, aeff_list, nhw_list, constructs a map from
    (A,Aeff) to the *.sh file that will be submitted to the cluster
    :param a_list: ordered list of mass numbers (A)
    :param aeff_list: ordered list of effective mass numbers (Aeff),
    corresponding to the mass numbers in a_list
    :param nhw_list: ordered list of Nhw corresponding to the mass numbers in
    a_list
    :param z: proton number (Z)
    :param n1: max allowed one-particle state (for TBME truncation)
    :param n2: max allowed two-particle state (for TBME truncation)
    :param a_aeff_to_dirpath_map: map from (A,Aeff) to the directory in which
    an NCSD calculation will be done
    :param scalefactor: value by which off-diagonal coupling terms in the
    TBME interaction were scaled.
    :param remove_protons: either true or false, indicating whether the
    proton parts of the interaction (Vpp and Vpn) are being scaled to 0
    """
    a_aeff_jobsub_map = dict()
    fname_fmt = FNAME_FMT_JOBSUB
    if remove_protons:
        fname_fmt += FNAME_ADD_JOBSUB_IPROT
    if scalefactor is not None:
        fname_fmt += FNAME_ADD_JOBSUB_SF
    fname_fmt += FNAME_EXT_JOBSUB
    for a, aeff, nhw in zip(a_list, aeff_list, nhw_list):
        args = (_get_name(z), a, aeff, nhw, n1, n2)
        if scalefactor is not None:
            args += (scalefactor,)
        a_aeff_jobsub_map[(a, aeff)] = path.join(
            a_aeff_to_dirpath_map[(a, aeff)], fname_fmt % args)
    a_aeff_to_ncsd = dict()
    a_aeff_to_vce = dict()
    for k, v in a_aeff_jobsub_map.items():
        a_aeff_to_ncsd[k] = v + '_NCSD'
        a_aeff_to_vce[k] = v + '_TRDENS'
    return a_aeff_to_ncsd, a_aeff_to_vce


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
        remove_protons, beta_cm, num_states, num_iter, dpath_templates,
        dpath_results, cluster_submit=False, walltime=None, progress=False,
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
    :param nshell: shell (0: s, 1: p, 2: sd, ...)
    :param n1: max allowed 1-particle state
    :param n2: max allowed 2-particle state
    :param scalefactor: factor by which to scale off-diagonal coupling terms of
    the TBME interaction
    :param remove_protons: if true, scales the proton-proton and
    proton-neutron parts of the interaction to 0
    :param beta_cm: ceter of mass beta term
    :param num_states: number of NCSD states to calculate
    :param num_iter: number of iterations for lanczos algorithm
    :param dpath_templates: path to the templates directory
    :param dpath_results: path to the results directory
    :param progress: if true, shows a progress bar
    :param cluster_submit: if true, submits the job to the OpenMP cluster
    using qsub
    :param walltime: walltime to be allotted to a cluster submission
    :param force: if true, forces re-truncation of the TBME interaction file
    :return (A,Aeff)->dir, (A,Aeff)->*.egv, (A,Aeff)->*.sh, (A,Aeff)->*.out
    """
    # get maps
    a_aeff_to_dir_map = get_a_aeff_to_dpath_map(
        a_list=a_list, aeff_list=aeff_list, nhw_list=nhw_list, z=z, n1=n1,
        n2=n2, scalefactor=scalefactor, remove_protons=remove_protons,
        dpath_results=dpath_results,
    )
    a_aeff_to_outfile_map = get_a_aeff_to_outfile_fpath_map(
        a_list=a_list, aeff_list=aeff_list, nhw_list=nhw_list, z=z, n1=n1,
        n2=n2, a_aeff_to_dirpath_map=a_aeff_to_dir_map, scalefactor=scalefactor,
        remove_protons=remove_protons,
    )
    a_aeff_to_egvfile_map = _get_a_aeff_to_egv_fpath_map(
        a_list=a_list, aeff_list=aeff_list, nhw_list=nhw_list,
        a_aeff_to_dirpath_map=a_aeff_to_dir_map
    )
    if cluster_submit:
        a_aeff_to_ncsd_map, a_aeff_to_vce_map = _get_a_aeff_to_jobsub_fpath_map(
            a_list=a_list, aeff_list=aeff_list, nhw_list=nhw_list, z=z, n1=n1,
            n2=n2, a_aeff_to_dirpath_map=a_aeff_to_dir_map,
            scalefactor=scalefactor, remove_protons=remove_protons,
        )
    else:
        a_aeff_to_ncsd_map = dict()
        a_aeff_to_vce_map = dict()
    # make stuff
    if progress:
        print '  Making directories'
    _make_base_directories(a_aeff_to_dpath_map=a_aeff_to_dir_map)
    if progress and (n1, n2) != (15, 15):
        print '  Truncating interaction to N1=%d N2=%d' % (n1, n2)
    if progress and (scalefactor is not None) and (scalefactor != 1.0):
        print '  Scaling interaction by %3.2f' % scalefactor
    fname_tbme, lname_tbme = _truncate_spaces(
        nshell=nshell, n1=n1, n2=n2, dirpaths=a_aeff_to_dir_map.values(),
        scalefactor=scalefactor, remove_protons=remove_protons, force=force
    )
    if progress:
        print '  Writing mfdp files'
    _make_mfdp_files(
        a_list=a_list, aeff_list=aeff_list, nhw_list=nhw_list,
        a_aeff_to_dpath_map=a_aeff_to_dir_map,
        a_aeff_to_outfile_fpath_map=a_aeff_to_outfile_map,
        fname_tbme=lname_tbme, dpath_temp=dpath_templates, z=z, nshell=nshell,
        n_1=n1, n_2=n2, beta_cm=beta_cm, num_states=num_states,
        num_iter=num_iter,
    )
    if cluster_submit:
        if progress:
            print '  Writing cluster submit files'
        _make_job_submit_files(
            a_aeff_to_jobsub_fpath_map=a_aeff_to_ncsd_map, walltime=walltime,
            fname_tmp_jobsub=FNAME_TMP_JOBSUB_NCSD
        )
        _make_job_submit_files(
            a_aeff_to_jobsub_fpath_map=a_aeff_to_vce_map, walltime=walltime,
            fname_tmp_jobsub=FNAME_TMP_JOBSUB_TRDENS
        )
    return (a_aeff_to_dir_map, a_aeff_to_egvfile_map, a_aeff_to_ncsd_map,
            a_aeff_to_vce_map, a_aeff_to_outfile_map)
