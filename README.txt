Replace your pipeline module with this file:

  pipeline/modules/spechla.nf  <= modules/spechla.nf from this zip

Key fix: escape bash command substitutions like \$(...) so Nextflow/Groovy compiles.
Also adds robust SpecHLA.sh discovery and fails loudly if results are missing.

Optional params you can add (YAML or CLI):
  spechla_home:   /path/to/SpecHLA
  spechla_script: /path/to/SpecHLA/script/whole/SpecHLA.sh
