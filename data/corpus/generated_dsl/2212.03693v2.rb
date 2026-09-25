# ============================================================================
#  Search for e+e- -> Omega- Omega+ with a single-hyperon tag
#  Eight c.m. energy points: 3.4900, 3.5080, 3.5097, 3.5104, 3.5146,
#                            3.5815, 3.6500, 3.6702 GeV
#
#  Two charge-conjugate tag selections (one Algorithm each):
#    I  Omega- tag : Omega- -> Lambda K-,  Lambda -> p pi-              (recoil -> Omega+)
#    II Omega+ tag : Omega+ -> anti-Lambda K+, anti-Lambda -> pbar pi+  (recoil -> Omega-)
#
#  The untagged (recoil) hyperon is inferred from the recoil four-momentum:
#  NO 4C kinematic fit is performed.
# ============================================================================

### ------------------------------ Datasets ------------------------------ ###
data_points = [
  DatasetManager.real_data.find("703_3490"),   # 3.4900 GeV
  DatasetManager.real_data.find("703_3508"),   # 3.5080 GeV
  DatasetManager.real_data.find("703_3509"),   # 3.5097 GeV
  DatasetManager.real_data.find("703_3510"),   # 3.5104 GeV
  DatasetManager.real_data.find("703_3514"),   # 3.5146 GeV
  DatasetManager.real_data.find("704_3581"),   # 3.5815 GeV
  DatasetManager.real_data.find("709_3650"),   # 3.6500 GeV
  DatasetManager.real_data.find("704_3670"),   # 3.6702 GeV
]

incMC_points = [
  DatasetManager.inclusive_mc.find("703_3490"),
  DatasetManager.inclusive_mc.find("703_3508"),
  DatasetManager.inclusive_mc.find("703_3509"),
  DatasetManager.inclusive_mc.find("703_3510"),
  DatasetManager.inclusive_mc.find("703_3514"),
  DatasetManager.inclusive_mc.find("704_3581"),
  DatasetManager.inclusive_mc.find("709_3650"),
  DatasetManager.inclusive_mc.find("704_3670"),
]

### ---------------- Signal decay card (EvtGen syntax) ---------------- ###
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.0 Omega- anti-Omega+ PHSP;
  Enddecay

  Decay Omega-
  1.0 Lambda0 K- PHSP;
  Enddecay

  Decay anti-Omega+
  1.0 anti-Lambda0 K+ PHSP;
  Enddecay

  Decay Lambda0
  1.0 p+ pi- PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0 anti-p- pi+ PHSP;
  Enddecay

  End
DECAYCARD

### ------- Exclusive MC: e+e- -> Omega- Omega+ (200k per energy point) ------- ###
# One signal sample is generated per scan point (same decay card, same size);
# both tag selections run on the full set of generated samples.
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_omega_omega"
  config.events        = 200000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### =================== Algorithm I : Omega- tag =================== ###
alg_name_minus = "OmegaMinusTag"
alg_minus = Algorithm.new(alg_name_minus)
alg_minus
  .set_header(["#{alg_name_minus}Alg/#{alg_name_minus}.h"])
  .set_constant({ "ECMS" => [:double, 3.510] })   # nominal c.m. energy fixed to 3.510 GeV

sel_minus = Selection.new
  .select_track {
      cos_theta 0.93      # |cos(theta)| < 0.93
      Vz        10.0      # |Vz| < 10 cm
      Vr        1.0       # Vr < 1 cm
      nChrp     ">=1"     # >= 1 positive track  (p+)
      nChrn     ">=2"     # >= 2 negative tracks (pi-, K-)
  }
  .pid(method: :probability) {          # probability PID from dE/dx and TOF
      prob_cut 0.001
      identify :pion,   against: [:kaon, :proton]    # pi+ / pi-
      identify :kaon,   against: [:pion, :proton]    # K+  / K-
      identify :proton, against: [:kaon, :pion]      # p   / pbar
      npim ">=1"
      nkm  ">=1"
      nprp ">=1"
  }
  .secondary_vertex_fit([:prp, :pim]) {   # Lambda -> p pi-
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .partial_miss([2]) {                    # miss anti-Omega+ (recID 2); reconstruct Omega- side
      best_combination_by_mass :Omega_minus, 1.67245   # Omega- candidate closest to 1.67245 GeV
      require_recoil_mass 1.59, 1.75                   # loose recoil-mass (Omega+) window
  }

# BOSS-side procedures with no dedicated DSL construct
alg_minus.note(:lambda_mass_and_decay_length,
               "Lambda candidate required to satisfy |M(p pi-) - m_Lambda| < 4 MeV and a " \
               "positive decay-length significance; the secondary-vertex-fit construct has " \
               "no mass-window or decay-length selector.")

alg_minus.with_decay_card(decay_card_signal).apply(sel_minus)

### ================ Algorithm II : Omega+ tag (c.c.) ================ ###
alg_name_plus = "OmegaPlusTag"
alg_plus = Algorithm.new(alg_name_plus)
alg_plus
  .set_header(["#{alg_name_plus}Alg/#{alg_name_plus}.h"])
  .set_constant({ "ECMS" => [:double, 3.510] })

sel_plus = Selection.new
  .select_track {
      cos_theta 0.93      # |cos(theta)| < 0.93
      Vz        10.0      # |Vz| < 10 cm
      Vr        1.0       # Vr < 1 cm
      nChrp     ">=2"     # >= 2 positive tracks (pi+, K+)
      nChrn     ">=1"     # >= 1 negative track  (pbar)
  }
  .pid(method: :probability) {
      prob_cut 0.001
      identify :pion,   against: [:kaon, :proton]
      identify :kaon,   against: [:pion, :proton]
      identify :proton, against: [:kaon, :pion]
      npip ">=1"
      nkp  ">=1"
      nprm ">=1"
  }
  .secondary_vertex_fit([:prm, :pip]) {   # anti-Lambda -> pbar pi+
      build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .partial_miss([1]) {                    # miss Omega- (recID 1); reconstruct anti-Omega+ side
      best_combination_by_mass :Omega_bar, 1.67245     # anti-Omega+ candidate closest to 1.67245 GeV
      require_recoil_mass 1.59, 1.75
  }

alg_plus.note(:lambda_mass_and_decay_length,
              "anti-Lambda candidate required to satisfy |M(pbar pi+) - m_Lambda| < 4 MeV and a " \
              "positive decay-length significance; not expressible in the secondary-vertex-fit DSL.")

alg_plus.with_decay_card(decay_card_signal).apply(sel_plus)

### ------------------------------ Execution ------------------------------ ###
root_files_minus = alg_minus.execute_on(data_points + incMC_points + exMCs_signal)
root_files_plus  = alg_plus.execute_on(data_points + incMC_points + exMCs_signal)