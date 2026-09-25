# Paper: 2307.15894v2
# e+e- → Sigma+ Sigma- at 2.3960, 2.6454, 2.9000 GeV
# ConExc generator (mode 44: Sigma- anti-Sigma+)
# Two methods: single-tag at 2.3960, double-tag with missing pi0 at higher energies

### Dataset preparation ###
# 2.6454 GeV = combined 2.6444 + 2.6464 datasets

# At 2.3960 GeV: single-tag method
data_2396 = DatasetManager.load_real_data.find("713_Rscan_2396")

# At 2.6454 GeV: combined 2.6444 + 2.6464 (not in a single dataset, use both)
data_2644 = DatasetManager.load_real_data.find("713_Rscan_2644")
data_2646 = DatasetManager.load_real_data.find("713_Rscan_2646")

# At 2.9000 GeV: double-tag method
data_2900 = DatasetManager.load_real_data.find("713_Rscan_2900")

all_datasets = [data_2396, data_2644, data_2646, data_2900]
all_incMC = all_datasets.map { |ds| DatasetManager.load_inclusive_mc.find(ds.sample_name) rescue nil }.compact

# ConExc decay card: mode 44 = Sigma- anti-Sigma+
# Particle vpho is auto-injected per energy point by DSL
conexc_card = <<~DECAYCARD
    Decay vpho
    1 ConExc 44;
    Enddecay
    Decay vhdr
    1 Sigma- anti-Sigma+ PHSP;
    Enddecay
    Decay anti-Sigma-
    1.0 anti-p- pi0 PHSP;
    Enddecay
    Decay Sigma+
    1.0 p+ pi0 PHSP;
    Enddecay
    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

exMC_2396 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "ee_SigmaSigma_2396"
  config.related_dataset = data_2396
  config.events = 6000000
  config.decay_card = conexc_card
  config.cross_section = :default
end

exMC_26454 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "ee_SigmaSigma_26454"
  config.related_dataset = data_2644
  config.events = 6000000
  config.decay_card = conexc_card
  config.cross_section = :default
end

exMC_2900 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "ee_SigmaSigma_2900"
  config.related_dataset = data_2900
  config.events = 6000000
  config.decay_card = conexc_card
  config.cross_section = :default
end

### Event selection (BOSS) ###
# Two separate analysis algorithms for the two methods

# ===== Algorithm 1: Single-tag at 2.3960 GeV =====
# Reconstruct anti-Sigma- → anti-proton pi0
# At least 1 antiproton, at least 2 photons
# Uses dE/dx only (tracks cannot reach TOF at this energy)

alg_st = Algorithm.new("EEtoSigmaSigma_SingleTag")
alg_st.set_header(["EEtoSigmaSigma_STAlg/EEtoSigmaSigma_ST.h"])
       .set_alias({"std::vector<double>" => "Vdouble"})
# No ECMS — scan energy

sel_st = Selection.new
sel_st.select_track {
            cos_theta   0.93
            Vz   10.0
            Vr   1.0
            nChrp   ">=0"
            nChrn   ">=1"     # At least one good antiproton
          }
         .select_photon {
            angle_to_track   20.0    # >20 degrees from nearest track
            energyThreshold_b   0.025
            energyThreshold_e   0.050
            nGam   ">=2"
          }
         .pid(method: :probability) {
            prob_cut   0.001
            identify :proton, against: [:kaon, :pion]
            nprm   ">=1"    # At least one antiproton
            # Note: at 2.3960 GeV, only dE/dx used (TOF inaccessible)
         }
         .remove([:prp <= :chrgp])
         .remove([:prm <= :chrgn])
         .select_isolated_photon {
            angle_to_prm_track   20.0
            nGam   ">=2"
         }
         # Reconstruct pi0 from gamma gamma, constrain to pi0 mass
         .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0   ">=1"
         }
         # Combine anti-proton + pi0 → anti-Sigma-
         .build_virtual_particle(:anti_Sigma_minus, from: [:prm, :pi0])
         # In ROOT: anti-Sigma- selected via DeltaE and Mbc variables
         # DeltaE = E_total - E_beam, Mbc = sqrt(E_beam^2/c^4 - P_Sigmabar^2/c^2)
         # Cuts: 0.126 < M_gammagamma < 0.139, -0.013 < DeltaE < 0.005
         # Signal region: 1.185 < Mbc < 1.191 GeV/c^2

alg_st.with_decay_card(conexc_card).apply(sel_st)
alg_st.note(:single_tag_method, "Single-tag: reconstruct only anti-Sigma- → anti-p pi0; Sigma+ side inclusive")
alg_st.note(:dEdx_only_pid, "At 2.3960 GeV only dE/dx used for PID (TOF inaccessible for low-momentum tracks)")
alg_st.note(:selection_cuts, "M_gammagamma in [0.126, 0.139] GeV; DeltaE in [-0.013, 0.005] GeV")
alg_st.note(:signal_region, "Mbc in [1.185, 1.191] GeV/c^2 for angular analysis")
alg_st.note(:angular_analysis, "Joint angular distribution fit for alpha and sin(DeltaPhi); sin(DeltaPhi) only accessible (one hyperon reconstructed)")
alg_st.note(:form_factors, "|GE/GM| = sqrt(s(1-alpha)/(4M_Sigma^2(1+alpha)))")

# ===== Algorithm 2: Double-tag with missing pi0 at 2.6454 and 2.9000 GeV =====
# proton + antiproton + pi0 (reconstructed) + missing pi0
# 2C kinematic fit: constrain reconstructed pi0 mass, missing pi0 free

alg_dt = Algorithm.new("EEtoSigmaSigma_DoubleTag")
alg_dt.set_header(["EEtoSigmaSigma_DTAlg/EEtoSigmaSigma_DT.h"])
       .set_alias({"std::vector<double>" => "Vdouble"})
# No ECMS — scan energy

sel_dt = Selection.new
sel_dt.select_track {
            cos_theta   0.93
            Vz   10.0
            Vr   1.0
            nChrp   ">=1"
            nChrn   ">=1"
          }
         .select_photon {
            angle_to_track   20.0
            energyThreshold_b   0.025
            energyThreshold_e   0.050
            nGam   ">=2"
          }
         .pid(method: :probability) {
            prob_cut   0.001
            identify :proton, against: [:kaon, :pion]
            nprp   ">=1"
            nprm   ">=1"
          }
         .remove([:prp <= :chrgp])
         .remove([:prm <= :chrgn])
         .select_isolated_photon {
            angle_to_prp_track   20.0
            angle_to_prm_track   20.0
            nGam   ">=2"
         }
         # Reconstruct one pi0 from gamma gamma
         .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0   ">=1"
         }
         # 2C kinematic fit: e+e- → p anti-p pi0(reconstructed) pi0(missing)
         # The reconstructed pi0 mass is constrained to nominal; the other pi0 is missing
         .kinematic_fit([:prp, :prm, :pi0]) {
            nominal
            constrain_four_momentum
            miss_track_of(:pi0)   # Allow one missing pi0
            chi2_cut 15            # chi2(2C) < 15
            # In ROOT: loop over pi0 candidates, pick best chi2 combination
            # In ROOT: pair pi0 with proton or antiproton → min |M(p/anti-p pi0) - M_Sigma+|
            # Signal region: 1.175 < M_Sigma_tag < 1.200 GeV/c^2
            # Recoil mass M_Sigma_rec used as observable
         }

alg_dt.with_decay_card(conexc_card).apply(sel_dt)
alg_dt.note(:double_tag_method, "Reconstruct p + anti-p + pi0, treat other pi0 as missing in 2C fit")
alg_dt.note(:missing_pi0_fit, "2C kinematic fit: constrain pi0 mass, total 4-momentum; missing pi0 has free 3-momentum")
alg_dt.note(:sigma_tag, "Best Sigma candidate: min |M(p/anti-p + pi0) - M_Sigma+|; signal region M_Sigma_tag in [1.175, 1.200]")
alg_dt.note(:recoil_mass, "Sigma0 inferred from recoil mass spectrum against Sigma_tag; signal extracted via unbinned ML fit")
alg_dt.note(:angular_analysis, "Full angular distribution fit for alpha and DeltaPhi from 5D angular variables")
alg_dt.note(:cp_conservation, "CP conservation assumed: alpha_Sigma+ = -alpha_Sigma- = -0.980")

root_files_st = alg_st.execute_on([data_2396])
root_files_dt = alg_dt.execute_on([data_2644, data_2646, data_2900])