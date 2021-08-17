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
params.input_folder = null
params.output_folder = "skat_output"

log.info ""
log.info "-----------------------------------------------------------------------"
log.info "run_skat.nf"
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
    log.info "--input_folder              FOLDER                 Folder containing one input file per window from the create_skat_data.nf script"
    log.info ""
    log.info "Optional arguments:"
    log.info '--output_folder             FOLDER                 Output folder (default: skat_output)'
    log.info ""
    log.info "Flags:"
    log.info "--help                                             Display this message"
    log.info ""
    exit 1
}

assert (params.input_folder != null) : "please provide the --input_folder option"

input_files = Channel.fromPath( params.input_folder+'/*.Rdata' )

process skat {

  input:
  file f from input_files

  output:
  file '*pvalue' into skatpvalues

  shell:
  '''
  Rscript !{baseDir}/bin/SKAT.R --input_file=!{f}
  '''
}

process manhattan {

  input:
  file f from skatpvalues.collect()

  output:
  file '*.pdf' into wind_list

  shell:
  '''
  Rscript !{baseDir}/bin/manhattan_plot.R --manhattan_plot_function=!{baseDir}/bin/manhattan_plot_function.R
  '''
}
