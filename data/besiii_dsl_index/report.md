# Fact index build report

- papers: 862
- chunks: 66360
- facts: 219289
- annotation status: {'ok': 853, 'invalid': 8, 'partial': 1}
- papers with zero facts: 0
- DSL parse failures: 0

## Facts by kind

- `table_cell`: 141946
- `text_percent`: 20033
- `dsl_param`: 16812
- `text_value`: 8445
- `dsl_method`: 7247
- `note_fact`: 7018
- `dsl_mc`: 6942
- `dsl_decay`: 6162
- `dsl_dataset`: 4684

## Selectivity highlights (corpus_facts.md gates)

- boilerplate `dsl_param|pid.prob_cut|0.001` in 453 papers
- boilerplate `dsl_param|select_track.cos_theta|0.93` in 341 papers
- boilerplate `dsl_param|kinematic_fit.chi2_cut|200` in 319 papers
- boilerplate `dsl_param|select_photon.tdc_emc_start|0` in 285 papers
- boilerplate `dsl_param|select_photon.energyThreshold_b|0.025` in 276 papers
- boilerplate `dsl_param|select_photon.energyThreshold_e|0.050` in 269 papers
- boilerplate `dsl_param|select_photon.tdc_emc_end|14` in 253 papers
- boilerplate `dsl_param|select_track.Vz|10.0` in 250 papers
- boilerplate `dsl_param|select_track.Vr|1.0` in 244 papers
- boilerplate `dsl_param|select_track.nNet|==0` in 210 papers
- boilerplate `dsl_param|select_photon.angle_to_track|10.0` in 202 papers
- boilerplate `dsl_param|tag_side.charm|-1` in 171 papers
- boilerplate `dsl_param|selection.cos_theta|0.93` in 165 papers
- boilerplate `dsl_param|kalman_kinematic_fit.chi2_cut|25` in 150 papers
- boilerplate `dsl_param|selection.Vr|1.0` in 150 papers

- datasets reused (≥2 papers): 151
- decays reused (≥2 papers): 400
