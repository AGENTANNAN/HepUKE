# Paper: 2308.03361v1
# e+e- → Lambda Sigma0 + c.c. at 14 energy points from 2.3094 to 3.0800 GeV
# ConExc generator (mode 4: Lambda anti-Sigma0; c.c. mode 5: Sigma0 anti-Lambda)
# Two methods: indirect at 2.3094 (near threshold), single Lambda-tag at higher energies

### Dataset preparation ###
# 14 energy points
data_23094 = DatasetManager.load_real_data.find("713_Rscan_2309")
data_23864 = DatasetManager.load_real_data.find("713_Rscan_2386")
data_23960 = DatasetManager.load_real_data.find("713_Rscan_2396")
data_25000 = DatasetManager.load_real_data.find("713_Rscan_2500")
data_26444 = DatasetManager.load_real_data.find("713_Rscan_2644")
data_26464 = DatasetManager.load_real_data.find("713_Rscan_2646")
data_27000 = DatasetManager.load_real_data.find("713_Rscan_2700")
data_28000 = DatasetManager.load_real_data.find("713_Rscan_2800")
data_29000 = DatasetManager.load_real_data.find("713_Rscan_2900")
data_29500 = DatasetManager.load_real_data.find("713_Rscan_2950")
data_29810 = DatasetManager.load_real_data.find("713_Rscan_2981")
data_30000 = DatasetManager.load_real_data.find("713_Rscan_3000")
data_30200 = DatasetManager.load_real_data.find("713_Rscan_3020")
data_30800 = DatasetManager.load_real_data.find("713_Rscan_3080")

all_datasets = [data_23094, data_23864, data_23960, data_25000, data_26444, data_26464,
                data_27000, data_28000, data_29000, data_29500, data_29810,
                data_30000, data_30200, data_30800]
all_incMC = all_datasets.map { |ds| DatasetManager.load_inclusive_mc.find(ds.sample_name) rescue nil }.compact

# ConExc decay card: mode 4 = Lambda anti-Sigma0 (c.c. included in mode 5)
# Particle vpho is auto-injected per energy point by DSL
conexc_card = <<~DECAYCARD
    Decay vpho
    1 ConExc 4;
    Enddecay
    Decay vhdr
    1 Lambda anti-Sigma0 PHSP;
    Enddecay
    Decay anti-Sigma0
    1.0 gamma anti-Lambda PHSP;
    Enddecay
    Decay Lambda
    1.0 p+ pi- HypWK;
    Enddecay
    Decay anti-Lambda
    1.0 anti-p- pi+ HypWK;
    Enddecay
    End
DECAYCARD

exMC_scan = DatasetManager.create_exclusive_mc_for(all_datasets) do |config|
  config.sample_name = "ee_to_LambdaSigma0"
  config.events = 500000
  config.decay_card = conexc_card
  config.cross_section = :default
end

### Event selection (BOSS) ###

# ===== Method 1: Near-threshold indirect search at 2.3094 GeV =====
# This method is largely inexpressible in the BOSS DSL framework.
# The analysis uses: pi+ + pi- from primary vertex + secondary tracks from
# antiproton annihilation in beam pipe, signal extracted from pi+ momentum fit.
# We capture the track-level cuts and document the rest via .note().

alg_thr = Algorithm.new("EEtoLambdaSigma0_Threshold")
alg_thr.set_header(["EEtoLambdaSigma0_ThrAlg/EEtoLambdaSigma0_Thr.h"])
       .set_alias({"std::vector<double>" => "Vdouble"})
# No ECMS — scan energy

sel_thr = Selection.new
sel_thr.select_track {
            cos_theta   0.93
            Vz   10.0
            Vr   1.0
            nChrp   ">=1"
            nChrn   ">=1"
          }
         .pid(method: :probability) {
            prob_cut   0.001
            identify :pion, against: [:kaon, :proton]
            npip   "==1"
            npim   "==1"
            # At 2.3094 GeV only dE/dx used (TOF inaccessible)
         }

alg_thr.with_decay_card(conexc_card).apply(sel_thr)
alg_thr.note(:indirect_method, "Indirect search: identify signal via pi+ and pi- tracks from Lambda/Sigma0 decay + secondary tracks from anti-proton beam pipe annihilation")
alg_thr.note(:pion_vertex, "Vertex fit of two pion tracks; transverse distance to beam < 2 cm")
alg_thr.note(:pion_momentum, "pi- momentum in (0.08, 0.12) GeV/c to suppress background; signal extracted from pi+ momentum fit (MC-convolved Gaussian)")
alg_thr.note(:antiproton_secondary, "At least 2 additional tracks from common vertex, transverse distance 1-5 cm (beam pipe), indicating anti-proton annihilation")
alg_thr.note(:beam_background, "Beam-associated background estimated from 2.1250 GeV data (below threshold)")
alg_thr.note(:lambda_lambda_background, "Physical background from e+e- → Lambda anti-Lambda estimated from MC")
alg_thr.note(:dEdx_only, "At 2.3094 GeV only dE/dx used for PID (TOF inaccessible for low-momentum tracks)")
alg_thr.note(:cross_section, "sigma_B = N_obs / (L × epsilon × (1+delta) × B(Lambda→ppi) × B(Sigma0→gamma Lambda) × B(Lambda→ppi))")

# ===== Method 2: Single Lambda-tag at higher energies (12 points) =====
# Reconstruct Lambda → p pi- via secondary vertex fit
# Infer Sigma0 from recoil mass
# 2.3864 to 3.0800 GeV

alg_tag = Algorithm.new("EEtoLambdaSigma0_SingleTag")
alg_tag.set_header(["EEtoLambdaSigma0_TagAlg/EEtoLambdaSigma0_Tag.h"])
       .set_alias({"std::vector<double>" => "Vdouble"})

sel_tag = Selection.new
sel_tag.select_track {
            cos_theta   0.93
            Vz   30.0     # Looser Vz cut (30 cm) for Lambda daughters
            Vr   10.0     # Looser Vr cut (10 cm) for Lambda daughters
            nChrp   ">=1"
            nChrn   ">=1"
          }
         .pid(method: :probability) {
            prob_cut   0.001
            identify :proton, against: [:kaon, :pion]
            identify :pion, against: [:kaon]
            nprp   "==1"
            npim   "==1"
            # Only one p pi- pair kept
            # In ROOT: highest likelihood method: L(h) > L(other_hypotheses) for all h
          }
         # Secondary vertex fit for Lambda → p pi-
         .secondary_vertex_fit([:prp, :pim]) {
            # Lambda candidate: decay length > 2 × vertex resolution
            # chi2(primary_vtx) + chi2(secondary_vtx) < 50
            # In ROOT: Lambda mass window [1.11, 1.12] GeV/c^2
            build_virtual_particle(:Lambda).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
         }
         # Recoil mass technique: infer Sigma0 from M_recoil against Lambda
         # M_recoil^2 = E_Sigma0^2/c^4 - |p_e+e- - p_Lambda|^2/c^2
         # where E_Sigma0 = (s + m_Sigma0^2 - m_Lambda^2) / (2*sqrt(s))

alg_tag.with_decay_card(conexc_card).apply(sel_tag)
alg_tag.note(:single_lambda_tag, "Single Lambda-tag: reconstruct only Lambda → p pi-; Sigma0 inferred from recoil mass spectrum")
alg_tag.note(:lambda_selection, "Lambda: secondary vertex fit, decay length > 2*sigma_vtx, chi2_sum < 50, M_ppi in [1.11, 1.12] GeV/c^2")
alg_tag.note(:recoil_mass, "M_recoil = sqrt( E_Sigma0^2/c^4 - |p_e+e- - p_Lambda|^2/c^2 ), E_Sigma0 = (s + m_Sigma0^2 - m_Lambda^2)/(2*sqrt(s))")
alg_tag.note(:signal_extraction, "Signal yield from unbinned ML fit to M_recoil: MC-convolved Gaussian signal + Lambda_antiLambda MC + Sigma0_Sigma0 MC + Lambda_from_Sigma0 + sideband background")
alg_tag.note(:lambda_sidebands, "Sidebands: M_ppi in [1.095, 1.105] and [1.125, 1.135] GeV/c^2")
alg_tag.note(:dominant_backgrounds, "e+e- → Lambda anti-Lambda and e+e- → Sigma0 anti-Sigma0")
alg_tag.note(:cross_section, "sigma_B = N_obs / (L × epsilon × (1+delta) × B(Lambda→ppi))")

# Execute: indirect method on 2.3094 only; single-tag on remaining 12 points
higher_energy_datasets = [data_23864, data_23960, data_25000, data_26444, data_26464,
                          data_27000, data_28000, data_29000, data_29500, data_29810,
                          data_30000, data_30200, data_30800]
higher_incMC = higher_energy_datasets.map { |ds| DatasetManager.load_inclusive_mc.find(ds.sample_name) rescue nil }.compact

root_files_thr = alg_thr.execute_on([data_23094])
root_files_tag = alg_tag.execute_on(higher_energy_datasets + higher_incMC + exMC_scan)