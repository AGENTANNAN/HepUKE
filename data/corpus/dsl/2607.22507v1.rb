# Dataset preparation
# Analysis uses 56 c.m. energy points between 3.510 and 4.951 GeV (44 fb^-1 total),
# spanning psi(3770), psi(4040), psi(4160), Y(4230), Y(4360), psi(4415), Y(4500),
# Y(4660), Y(4710). Use psi(3770) sample as the reference; the same analysis is
# executed on every scan point (see note below).
data_ref  = DatasetManager.real_data.find("712_3773")
incMC_ref = DatasetManager.inclusive_mc.find("712_3773")

# Signal decay card: e+e- -> K_S0 anti-Xi+ Sigma- (charge conjugate implied).
# Reconstructed side: K_S0 -> pi+ pi-, anti-Xi+ -> anti-Lambda pi+ with anti-Lambda -> anti-p pi+.
# Sigma- is inferred from the recoil against the reconstructed K_S0 anti-Xi+ system
# (partial reconstruction).
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.0000 K_S0 anti-Xi+ Sigma-           PHSP;
  Enddecay

  Decay anti-Xi+
  1.0000 anti-Lambda0 pi+               PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+                    PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi-                        PHSP;
  Enddecay

  Decay Sigma-
  1.0000 n0 pi-                         PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "ee_KsXibarSigmaMinus"
  config.related_dataset = data_ref
  config.events          = 400000     # 4e5 PHSP events per energy point
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Event selection (BOSS)
alg_name = "eeKsXibarSigma"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 3.773] })

event_selection = Selection.new
event_selection.select_track {
                 cos_theta 0.93           # |cos theta| < 0.93 for MDC tracks
                 Vz    100.0
                 Vr    10.0
                 nChrp ">=3"              # at least 3 pi+ from Ks (1) + Lambda_bar (1) + Xi_bar (1)
                 nChrn ">=2"              # at least 2 negative tracks (anti-p + pi-)
               }
               .pid(method: :probability) {
                 prob_cut 0.001
                 # Proton hypothesis: L(p) > L(K) and L(p) > L(pi)
                 identify :proton, against: [:kaon, :pion]
                 # Pion hypothesis: L(pi) > L(K) and L(pi) > L(p)
                 identify :pion,   against: [:kaon, :proton]
                 nprm ">=1"               # at least 1 anti-proton
                 npip ">=3"               # at least 3 pi+
                 npim ">=1"               # at least 1 pi-
               }
               # Reconstruct K_S0 -> pi+ pi- via secondary vertex fit
               .secondary_vertex_fit([:pip, :pim]) {
                 build_virtual_particle(:K_S0).by_minimizing_mass_difference
                 remove_used_particle_from_candidate_list
               }
               # Reconstruct anti-Lambda -> anti-p pi+ via secondary vertex fit
               .secondary_vertex_fit([:prm, :pip]) {
                 build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                 remove_used_particle_from_candidate_list
               }
               # Reconstruct anti-Xi+ -> anti-Lambda pi+ via secondary vertex fit
               .secondary_vertex_fit([:Lambda_bar, :pip]) {
                 build_virtual_particle(:"anti-Xi+").by_minimizing_mass_difference
                 remove_used_particle_from_candidate_list
               }
               # Partial reconstruction: Sigma- inferred from recoil against K_S0 anti-Xi+ system
               .partial_miss([]) {
                 require_recoil_mass 1.1, 1.3
               }

alg.note(:ks_mass_window,
         "|M(pi+pi-) - m_KS0| < 12 MeV/c^2 required after the K_S0 secondary vertex fit.")
   .note(:ks_decay_length,
         "K_S0 decay length significance L/DeltaL > 2 required.")
   .note(:lambda_bar_mass_window,
         "|M(pbar pi+) - m_Lambda_bar| < 5 MeV/c^2 and decay length > 0 required.")
   .note(:xi_bar_mass_window,
         "|M(Lambda_bar pi+) - m_Xi_bar+| < 8 MeV/c^2 and decay length > 0 required.")
   .note(:best_combination,
         "Best anti-Lambda and anti-Xi+ combination chosen by minimising " \
         "delta_min = sqrt(|M(pbar pi+) - m_Lambda|^2 + |M(Lambda_bar pi+) - m_Xi|^2).")
   .note(:multi_energy_scan,
         "Analysis runs on 56 e+e- c.m. energy points between 3.510 and 4.951 GeV " \
         "(sum 44 fb^-1). The same selection is executed on every scan point; the reference " \
         "dataset here (712_3773) is a placeholder representative of the psi(3770) energy. " \
         "For a production scan, use DatasetManager.create_exclusive_mc_for(scan_points).")
   .note(:sigma_minus_recoil,
         "Sigma- reconstructed from the recoil against the K_S0 anti-Xi+ system: " \
         "M_recoil = sqrt((sqrt(s) - E_KsXibar)^2 - |p_KsXibar|^2); signal window in [1.1, 1.3] GeV/c^2.")

alg.with_decay_card(decay_card_signal).apply(event_selection)
alg.execute_on([data_ref, incMC_ref, exMC_signal])
