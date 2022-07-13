#!/usr/bin/env nextflow

// Copyright (C) 2020 IRB Barcelona

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.help = null
params.df_windows = null
params.somatic_folder = null
params.somatic_files = null
params.germline_VCF = null
params.output_folder = "input_data"
params.nwindow_list = 120
params.bed_overlap = null

log.info ""
log.info "-----------------------------------------------------------------------"
log.info "create_skat_data.nf"
log.info "-----------------------------------------------------------------------"
log.info "Copyright (C) IRB Barcelona"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "                     USAGE                              "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run main.nf [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--df_windows                R file                 File .Rdata containing the dataframe of the windows"
    log.info "--somatic_folder            FOLDER                 Folder containing the somatic files"
    log.info "--somatic_files             FILE                   Text file containing the list of somatic files to analyze"
    log.info "--germline_VCF              FILE                   VCF file containing the genotypes"
    log.info ""
    log.info "Optional arguments:"
    log.info '--output_folder             FOLDER                 Output folder (default: input_data)'
    log.info "--nwindow_list              INT                    Number of chunks of windows to be run in parallel (default=120)"
    log.info "--bed_overlap               FILE                   BED file containing positions to extract windows, in list, that are overlapping"
    log.info '--mem                       INT                    Memory to allocate'
    log.info ""
    log.info "Flags:"
    log.info "--help                                             Display this message"
    log.info ""
    exit 1
}

assert (params.df_windows != null) : "please provide the --df_windows option"
assert (params.somatic_folder != null) : "please provide the --somatic_folder option"
assert (params.somatic_files != null) : "please provide the --somatic_files option"
assert (params.germline_VCF != null) : "please provide the --germline_VCF option"

log.info "Somatic folder: ${params.somatic_folder}"
log.info "Germline VCF: ${params.germline_VCF}"
log.info "Output folder: ${params.output_folder}"
log.info "Rdata windows: ${params.df_windows}"

somatic_files = file(params.somatic_files)
somatic_folder = file(params.somatic_folder)
germline_VCF = file(params.germline_VCF)
germline_VCF_tbi = file(params.germline_VCF + ".tbi")
df_windows = file(params.df_windows)

process output_windows {

  input:
  file df_windows

  output:
  file 'window*' into wind_list mode flatten

  shell:
  if (params.bed_overlap=="null") { par="--nwindow_list=!{params.nwindow_list}" } else { par="--bed_overlap=!{params.bed_overlap}" }
  '''
  Rscript !{baseDir}/bin/output_windows.R --df_windows=!{df_windows} !{par} 
  '''
}

process load_somatic {

  input:
  file somatic_files
  file somatic_folder

  output:
  file 'gr_mut*' into somatic_all

  shell:
  '''
  Rscript !{baseDir}/bin/load_somatic.R --somatic_files=!{somatic_files} --somatic_folder=!{somatic_folder}
  '''
}

process create_input {

  publishDir params.output_folder, mode: 'copy', pattern: "*.Rdata"

  tag {tag}

  input:
  file wlist from wind_list
  file somatic_Rdata from somatic_all 
  file somatic_files
  file somatic_folder
  file germline_VCF
  file germline_VCF_tbi

  output:
  file '*.Rdata'  optional true into skat_input

  shell:
  tag = wlist
  '''
  Rscript !{baseDir}/bin/create_input_data.R --wlist=!{wlist} --somatic_Rdata=!{somatic_Rdata} --somatic_files=!{somatic_files} --somatic_folder=!{somatic_folder} --germline_VCF=!{germline_VCF}
  '''
}
