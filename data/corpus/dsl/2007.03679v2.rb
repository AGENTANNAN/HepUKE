# BOSS Ruby DSL for arXiv:2007.03679v2
# "Model independent determination of the spin of the Omega- and its polarization
#  alignment in psi(3686) -> Omega- antiOmega+"
#
# Analysis: psi(3686) -> Omega- antiOmega+
#   Single-tag method: reconstruct one Omega and infer the other from missing mass.
#   Two tagging strategies: Omega- tagged (Omega- -> K- Lambda -> K- p pi-)
#                       and Omega-bar+ tagged (antiOmega+ -> K+ antiLambda -> K+ anti-p pi+)
#   Recoil mass window [1.640, 1.692] GeV/c2 selects the other Omega.
#   Angular distribution fit performed in ROOT on selected events.

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card: psi(3686) -> Omega- antiOmega+
# MC generated with PHSP, then weighted with measured helicity amplitudes for efficiency
decay_card = <<~DECAYCARD
    Decay psi(2S)
    1.000 Omega- anti-Omega+ PHSP;
    Enddecay

    Decay Omega-
    1.000 K- Lambda0 PHSP;
    Enddecay

    Decay anti-Omega+
    1.000 K+ anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay
End
DECAYCARD

# Exclusive signal MC for efficiency determination
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "omega_spin_signal_mc"
  config.related_dataset = psip_data
  config.events = 1_000_000
  config.decay_card = decay_card
  config.cross_section = :default
end

# ===========================================================================
# Algorithm 1: Omega- tagged
#   Reconstruct Omega- -> K- Lambda -> K- p pi-
#   Miss anti-Omega+ via recoil mass
# ===========================================================================
alg_omega = Algorithm.new("OmegaSpinOmega")
alg_omega.set_header(["OmegaSpinOmegaAlg/OmegaSpinOmega.h"])
          .set_constant({"ECMS" => [:double, 3.686]})
          .note(:pid_highest_probability, "PID assigns each track the hypothesis with highest probability; DSL uses prob_cut threshold which is an approximation")
          .note(:lambda_mass_window, "M(ppi-) mass window [1.110, 1.122] GeV/c2 applied after Lambda vertex fit; DSL expresses best-by-mass only")
          .note(:omega_mass_window, "M(K-Lambda) mass window [1.663, 1.681] GeV/c2 applied after Omega vertex fit; DSL best_combination_by_mass approximates this")

sel_omega = Selection.new
sel_omega.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=1"
  nChrn ">=2"
end
# PID: identify proton first, then kaon from remaining tracks
.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  nprp "==1"
end
.remove([:prp <= :chrgp])
.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  nkm "==1"
end
.remove([:km <= :chrgn])
# Remaining negative tracks assumed to be pions
.assign({chrgn: :pim})
# Transverse momentum cuts: p_t > 0.2 (p), 0.1 (K-), 0.05 (pi-) GeV/c
.for_each(:prp) do
  where { px.__pow__(2) + py.__pow__(2) < 0.04 }
  remove
end
.for_each(:km) do
  where { px.__pow__(2) + py.__pow__(2) < 0.01 }
  remove
end
.for_each(:pim) do
  where { px.__pow__(2) + py.__pow__(2) < 0.0025 }
  remove
end
# Lambda -> p pi- secondary vertex fit
.secondary_vertex_fit([:prp, :pim]) do
  build_virtual_particle(:Lambda).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end
# Partial reconstruction: miss anti-Omega+ branch; reconstruct Omega- from K- + Lambda
.partial_miss([2]) do
  best_combination_by_mass :Omega, 1.67245
  require_recoil_mass 1.640, 1.692
end

alg_omega.with_decay_card(decay_card).apply(sel_omega)

# ===========================================================================
# Algorithm 2: Omega-bar+ tagged
#   Reconstruct anti-Omega+ -> K+ anti-Lambda -> K+ anti-p pi+
#   Miss Omega- via recoil mass
# ===========================================================================
alg_omega_bar = Algorithm.new("OmegaSpinOmegaBar")
alg_omega_bar.set_header(["OmegaSpinOmegaBarAlg/OmegaSpinOmegaBar.h"])
             .set_constant({"ECMS" => [:double, 3.686]})
             .note(:pid_highest_probability, "PID assigns each track the hypothesis with highest probability; DSL uses prob_cut threshold which is an approximation")
             .note(:lambda_bar_mass_window, "M(anti-p pi+) mass window [1.110, 1.122] GeV/c2 applied after anti-Lambda vertex fit; DSL expresses best-by-mass only")
             .note(:omega_bar_mass_window, "M(K+ anti-Lambda) mass window [1.663, 1.681] GeV/c2 applied after anti-Omega vertex fit; DSL best_combination_by_mass approximates this")

sel_omega_bar = Selection.new
sel_omega_bar.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"
  nChrn ">=1"
end
# PID: identify anti-proton first, then anti-kaon from remaining tracks
.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  nprm "==1"
end
.remove([:prm <= :chrgn])
.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  nkp "==1"
end
.remove([:kp <= :chrgp])
# Remaining positive tracks assumed to be pions
.assign({chrgp: :pip})
# Transverse momentum cuts: p_t > 0.2 (anti-p), 0.1 (K+), 0.05 (pi+) GeV/c
.for_each(:prm) do
  where { px.__pow__(2) + py.__pow__(2) < 0.04 }
  remove
end
.for_each(:kp) do
  where { px.__pow__(2) + py.__pow__(2) < 0.01 }
  remove
end
.for_each(:pip) do
  where { px.__pow__(2) + py.__pow__(2) < 0.0025 }
  remove
end
# anti-Lambda -> anti-p pi+ secondary vertex fit
.secondary_vertex_fit([:prm, :pip]) do
  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end
# Partial reconstruction: miss Omega- branch; reconstruct anti-Omega+ from K+ + anti-Lambda
.partial_miss([1]) do
  best_combination_by_mass :Omega_bar, 1.67245
  require_recoil_mass 1.640, 1.692
end

alg_omega_bar.with_decay_card(decay_card).apply(sel_omega_bar)

# Execute both algorithms on data and MC
alg_omega.execute_on([psip_data, psip_incMC, exMC_signal])
alg_omega_bar.execute_on([psip_data, psip_incMC, exMC_signal])