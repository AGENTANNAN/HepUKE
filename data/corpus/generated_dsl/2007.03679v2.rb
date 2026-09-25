# ============================================================
# psi(3686) -> Omega- anti-Omega+  single-tag analysis (BOSS part)
# sqrt(s) = 3.686 GeV; reconstruct one Omega and infer the recoiling
# Omega from the missing (recoil) mass.
# ============================================================

### Dataset description ###
data_3686  = DatasetManager.real_data.find("709_3686")      # psi(3686) real data at 3.686 GeV
incMC_3686 = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC

# Decay card: psi(3686) -> Omega- anti-Omega+, Omega- -> K- Lambda0,
# Lambda0 -> p+ pi-  (and the charge-conjugate chain).
# Signal is generated initially phase-space (PHSP); the measured
# helicity-amplitude reweighting is applied afterwards for efficiency.
decay_card_OmegaOmega = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Omega- anti-Omega+ PHSP;
    Enddecay

    Decay Omega-
    1.0000 K- Lambda0 PHSP;
    Enddecay

    Decay anti-Omega+
    1.0000 K+ anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC: 1,000,000 events of psi(3686) -> Omega- anti-Omega+
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_OmegaOmega"
  config.related_dataset = data_3686
  config.events          = 1_000_000
  config.decay_card      = decay_card_OmegaOmega
  config.cross_section   = :default
end

### ================= Event selection: Omega- single tag ================= ###
alg_name_om = "OmegaMinusTag"
alg_om = Algorithm.new(alg_name_om)
alg_om.set_header(["#{alg_name_om}Alg/#{alg_name_om}.h"])
      .set_constant({"ECMS" => [:double, 3.686]})
      .note(:signal_mc_reweighting,
            "exclusive signal MC generated phase-space (PHSP) and then reweighted by " \
            "the measured helicity amplitudes to obtain the efficiency")

sel_om = Selection.new
sel_om.select_track {
          cos_theta 0.93          # |cos(theta)| < 0.93
          Vz        10.0          # |Vz| < 10 cm
          Vr        1.0           # Vr < 1 cm
          nChrp     ">=1"         # at least one positive track
          nChrn     ">=2"         # at least two negative tracks
        }
      .pid(method: :probability) {
          prob_cut 0.001                                        # probability threshold 0.001
          identify :proton, against: [:kaon, :pion]             # one proton (p+) vs K/pi
          identify :kaon,   against: [:pion, :proton]           # one K- vs pi/p
          nprp "==1"
          nkm  "==1"
        }
      .remove([:prp <= :chrgp])     # remove identified proton from the positive list
      .remove([:km  <= :chrgn])     # remove identified K- from the negative list
      .assign({:chrgn => :pim})     # remaining negative track is taken as pi-
      .for_each(:prp) do
        where { pt_of(:prp) < 0.2 }     # pT(proton) > 0.2 GeV
        remove
      end
      .for_each(:km) do
        where { pt_of(:km) < 0.1 }      # pT(K-) > 0.1 GeV
        remove
      end
      .for_each(:pim) do
        where { pt_of(:pim) < 0.05 }    # pT(pi-) > 0.05 GeV
        remove
      end
      .secondary_vertex_fit([:prp, :pim]) {           # Lambda -> p pi- secondary vertex
          build_virtual_particle(:Lambda).by_minimizing_mass_difference
          remove_used_particle_from_candidate_list
      }
      # Partial reconstruction: tag Omega- (K- Lambda); the recoiling anti-Omega+
      # is inferred from the missing mass. No global nC kinematic fit / chi2 cut.
      .partial_rec([1, 3, 4, 5, 6]) do
          best_combination_by_mass :Omega_minus, 1.67245   # K-Lambda closest to nominal Omega mass
          require_recoil_mass 1.640, 1.692                 # recoil anti-Omega mass window
      end

alg_om.with_decay_card(decay_card_OmegaOmega).apply(sel_om)
alg_om.execute_on([data_3686, incMC_3686, exMC_signal])

### ================ Event selection: anti-Omega+ single tag (charge conjugate) ================ ###
alg_name_op = "OmegaPlusTag"
alg_op = Algorithm.new(alg_name_op)
alg_op.set_header(["#{alg_name_op}Alg/#{alg_name_op}.h"])
      .set_constant({"ECMS" => [:double, 3.686]})
      .note(:signal_mc_reweighting,
            "exclusive signal MC generated phase-space (PHSP) and then reweighted by " \
            "the measured helicity amplitudes to obtain the efficiency")

sel_op = Selection.new
sel_op.select_track {
          cos_theta 0.93          # |cos(theta)| < 0.93
          Vz        10.0          # |Vz| < 10 cm
          Vr        1.0           # Vr < 1 cm
          nChrp     ">=2"         # at least two positive tracks
          nChrn     ">=1"         # at least one negative track
        }
      .pid(method: :probability) {
          prob_cut 0.001                                        # probability threshold 0.001
          identify :proton, against: [:kaon, :pion]             # one anti-proton (pbar) vs K/pi
          identify :kaon,   against: [:pion, :proton]           # one K+ vs pi/p
          nprm "==1"
          nkp  "==1"
        }
      .remove([:prm <= :chrgn])     # remove identified anti-proton from the negative list
      .remove([:kp  <= :chrgp])     # remove identified K+ from the positive list
      .assign({:chrgp => :pip})     # remaining positive track is taken as pi+
      .for_each(:prm) do
        where { pt_of(:prm) < 0.2 }     # pT(anti-proton) > 0.2 GeV
        remove
      end
      .for_each(:kp) do
        where { pt_of(:kp) < 0.1 }      # pT(K+) > 0.1 GeV
        remove
      end
      .for_each(:pip) do
        where { pt_of(:pip) < 0.05 }    # pT(pi+) > 0.05 GeV
        remove
      end
      .secondary_vertex_fit([:prm, :pip]) {           # anti-Lambda -> pbar pi+ secondary vertex
          build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
          remove_used_particle_from_candidate_list
      }
      # Partial reconstruction: tag anti-Omega+ (K+ anti-Lambda); the recoiling
      # Omega- is inferred from the missing mass. No global nC kinematic fit / chi2 cut.
      .partial_rec([2, 7, 8, 9, 10]) do
          best_combination_by_mass :Omega_plus, 1.67245    # K+anti-Lambda closest to nominal Omega mass
          require_recoil_mass 1.640, 1.692                 # recoil Omega mass window
      end

alg_op.with_decay_card(decay_card_OmegaOmega).apply(sel_op)
alg_op.execute_on([data_3686, incMC_3686, exMC_signal])