# Paper: 2206.10900v2
# Title: Observation of Xi- transverse polarization in psi(3686) -> Xi- Xibar+
# Energy: 3.686 GeV (psi(3686), single energy point)
# Final state: Xi- -> Lambda pi-, Lambda -> p pi-; Xibar+ -> Lambdabar pi+, Lambdabar -> pbar pi+
# Full hyperon decay chain: 6 charged tracks: p, pbar, pi-, pi-, pi+, pi+

### Dataset preparation ###
data_709_3686 = DatasetManager.load_real_data.find("709_3686")
incMC_709_3686 = DatasetManager.load_inclusive_mc.find("709_3686")

all_data = [data_709_3686]
all_incMC = [incMC_709_3686]

# Decay card: psi(3686) -> Xi- Xibar+
# Xi- -> Lambda pi-, Lambda -> p pi-
# Xibar+ -> Lambdabar pi+, Lambdabar -> pbar pi+
decay_card = <<~DECAYCARD
    Decay psi(3686)
    1.0000  Xi-  anti-Xi-  anti-Xi+    HELAMP 1.0 0.0 0.693 0.0;
    Enddecay

    Decay Xi-
    1.0000  Lambda0  pi-               PHSP;
    Enddecay

    Decay anti-Xi-
    1.0000  anti-Lambda0  pi+          PHSP;
    Enddecay

    Decay anti-Xi+
    1.0000  anti-Lambda0  pi+          PHSP;
    Enddecay

    Decay Lambda0
    1.0000  p+  pi-                    PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000  anti-p-  pi+               PHSP;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "psi3686_XiXibar_signal"
  config.events         = 500000
  config.decay_card     = decay_card
  config.cross_section  = :default
end

### Event selection ###
alg = Algorithm.new("XiXibar_Polarization")
alg.set_header(["XiXibarPolAlg/XiXibarPol.h"])
   .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new

event_selection.select_track {
                 cos_theta 0.93
                 nTot ">=6"       # p, pbar, pi-, pi-, pi+, pi+
               }
               # PID: track momentum-based assignment
               # p > 0.5 GeV/c -> proton; p < 0.5 GeV/c -> pion
               .pid(method: :probability) {
                 prob_cut 0.001
                 identify :proton, against: [:kaon, :pion]
                 nprp ">=1"
                 nprm ">=1"
               }
               .remove([:prp <= :chrgp])
               .remove([:prm <= :chrgn])
               # Assign remaining charged tracks as pions
               .assign({:chrgp => :pip, :chrgn => :pim})
               # Require at least 2 pi+ and 2 pi-
               # Build Lambda candidates via secondary vertex fit of p pi-
               .secondary_vertex_fit([:prp, :pim]) {
                 build_virtual_particle(:Lambda).by_minimizing_mass_difference
                 remove_used_particle_from_candidate_list
               }
               # Build anti-Lambda candidates via secondary vertex fit of pbar pi+
               .secondary_vertex_fit([:prm, :pip]) {
                 build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                 remove_used_particle_from_candidate_list
               }
               # Build Xi- candidates via secondary vertex fit of Lambda pi-
               .secondary_vertex_fit([:Lambda, :pim]) {
                 build_virtual_particle(:Xi_minus).by_minimizing_mass_difference
                 remove_used_particle_from_candidate_list
               }
               # Build anti-Xi+ candidates via secondary vertex fit of Lambdabar pi+
               .secondary_vertex_fit([:Lambda_bar, :pip]) {
                 build_virtual_particle(:Xi_plus).by_minimizing_mass_difference
                 remove_used_particle_from_candidate_list
               }
               # 4C kinematic fit: four-momentum constraint to CMS
               .kinematic_fit([:Xi_minus, :Xi_plus]) {
                 nominal
                 constrain_four_momentum
                 chi2_cut 200
               }

alg.with_decay_card(decay_card).apply(event_selection)

# Note: the paper uses best-combination selection by minimizing mass differences
alg.note(:best_combination,
  "Best Lambda candidate by minimizing sqrt((M_pi- - m_Lambda)^2 + (M_pbar_pi+ - m_Lambda)^2). Best Xi candidate by minimizing sqrt((M_Lambda_pi- - m_Xi)^2 + (M_Lambdabar_pi+ - m_Xi)^2). Applied in ROOT.")

# Lambda mass window: 5 MeV/c^2 around nominal mass, optimized by S/sqrt(S+B) FOM
alg.note(:lambda_mass_window,
  "M(p pi-) within 5 MeV/c^2 of Lambda mass; optimized by FOM. Decay length > 0. Applied in ROOT.")

# Xi mass window: 8 MeV/c^2 around nominal Xi mass, optimized by FOM
alg.note(:xi_mass_window,
  "M(Lambda pi-) within 8 MeV/c^2 of Xi mass; optimized by FOM. Decay length > 0. Applied in ROOT.")

# Chi2_4C optimized by FOM: chi2_4C < 200
alg.note(:chi2_4C_optimization,
  "chi2_4C < 200 retained based on FOM optimization. Applied in ROOT.")

# Background from sideband regions; negligible (~0.6%)
alg.note(:background,
  "Background ~0.6% from non-Xi events (e+e- -> pi+pi-LambdaLambdabar) and continuum. Sideband method used. Applied in ROOT.")

# Angular distribution fit for polarization parameters is ROOT-only
alg.note(:angular_fit,
  "Unbinned maximum likelihood fit for Xi polarization parameters (alpha_psi, DeltaPhi, alpha_Xi, phi_Xi) from joint angular distribution. Applied in ROOT.")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on(all_datasets)