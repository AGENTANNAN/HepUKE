### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")       # psi(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # matching inclusive MC sample

# Decay card: psi(3686) -> K- Lambda anti-Xi+  (non-radiative mode)
decay_card_KLamXiBar = <<~DECAYCARD
    Decay psi(2S)
    1.0000 K- Lambda0 anti-Xi+         PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                      HypWK;
    Enddecay

    Decay anti-Xi+
    1.0000 anti-Lambda0 pi+            HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+                 HypWK;
    Enddecay

    End
DECAYCARD

# Decay card: psi(3686) -> gamma K- Lambda anti-Xi+  (radiative mode)
decay_card_gamKLamXiBar = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma K- Lambda0 anti-Xi+   PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                      HypWK;
    Enddecay

    Decay anti-Xi+
    1.0000 anti-Lambda0 pi+            HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+                 HypWK;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 200k events for each of the two decay modes
exMC_KLamXiBar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_KLamXiBar"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_KLamXiBar
  config.cross_section   = :default
end

exMC_gamKLamXiBar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamKLamXiBar"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_gamKLamXiBar
  config.cross_section   = :default
end

### Event selection (BOSS) - Mode I: psi(3686) -> K- Lambda anti-Xi+ ###
alg_name_KLamXiBar = "PsiPToKLamXiBar"
alg_KLamXiBar = Algorithm.new(alg_name_KLamXiBar)
alg_KLamXiBar.set_header(["#{alg_name_KLamXiBar}Alg/#{alg_name_KLamXiBar}.h"])
             .set_constant({"ECMS" => [:double, 3.686]})   # ECMS = 3.686 GeV
             .note(:kaon_highest_pid_cl,
                   "Only the K- candidate with the highest PID confidence level is kept as the "
                   "bachelor kaon; additional kaon candidates are treated as pions and additional "
                   "protons are likewise treated as pions. Approximated in the generated code by "
                   "retaining all kaon candidates in the kaon lists while the remaining tracks are "
                   "assigned to the pion lists, the best combination being chosen by the "
                   "kinematic fit.")
             .note(:kaon_from_ip,
                   "The bachelor K- is required to originate from the interaction point (IP); "
                   "this impact-parameter / primary-vertex requirement on the K- has no dedicated "
                   "DSL construct and is applied in the generated BOSS code.")

sel_KLamXiBar = Selection.new
  .select_track {                       # Charged track selection
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz        10.0                      # |Vz| < 10 cm  (= 100 mm)
    Vr        1.0                       # Vr < 1 cm     (= 10 mm)
    nChrp     ">=3"                     # at least 3 positively charged tracks
    nChrn     ">=3"                     # at least 3 negatively charged tracks
    nNet      "==0"                     # net charge zero
  }
  .pid(method: :probability) {          # PID by the probability method
    prob_cut 0.001                      # CL > 0.001
    identify :proton, against: [:kaon, :pion]   # pi/K/p separation (p+ and anti-p-)
    identify :kaon,   against: [:pion, :proton] # pi/K/p separation (K+ and K-)
    nprp ">=1"                          # at least one proton
    nprm ">=1"                          # at least one anti-proton
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])     # remove identified (anti-)protons from charge lists
  .assign({:chrgp => :pip, :chrgn => :pim})     # remaining tracks treated as pions
  # Secondary vertex fits (mass-difference minimisation) for Lambda, anti-Lambda and anti-Xi+
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:Lambda_bar, :pip]) {
    build_virtual_particle(:Xi_bar).by_minimizing_mass_difference   # anti-Xi+ -> anti-Lambda pi+
    remove_used_particle_from_candidate_list
  }
  # 4C kinematic fit: chi2 < 200 and Lambda / anti-Lambda mass windows
  .kinematic_fit([:km, :Lambda, :Lambda_bar, :Xi_bar]) {
    nominal
    constrain_four_momentum             # 4C energy-momentum constraint
    chi2_cut 200
    invariant_mass_of(:Lambda).within(1.110, 1.121)       # Lambda mass window (GeV/c^2)
    invariant_mass_of(:Lambda_bar).within(1.110, 1.121)   # anti-Lambda mass window (GeV/c^2)
  }

alg_KLamXiBar.with_decay_card(decay_card_KLamXiBar).apply(sel_KLamXiBar)

### Event selection (BOSS) - Mode II: psi(3686) -> gamma K- Lambda anti-Xi+ ###
alg_name_gamKLamXiBar = "PsiPToGamKLamXiBar"
alg_gamKLamXiBar = Algorithm.new(alg_name_gamKLamXiBar)
alg_gamKLamXiBar.set_header(["#{alg_name_gamKLamXiBar}Alg/#{alg_name_gamKLamXiBar}.h"])
                .set_constant({"ECMS" => [:double, 3.686]})
                .note(:kaon_highest_pid_cl,
                      "Only the K- candidate with the highest PID confidence level is kept as the "
                      "bachelor kaon; additional kaon candidates are treated as pions and additional "
                      "protons are likewise treated as pions. Approximated in the generated code by "
                      "retaining all kaon candidates in the kaon lists while the remaining tracks are "
                      "assigned to the pion lists, the best combination being chosen by the "
                      "kinematic fit.")
                .note(:kaon_from_ip,
                      "The bachelor K- is required to originate from the interaction point (IP); "
                      "this impact-parameter / primary-vertex requirement on the K- has no dedicated "
                      "DSL construct and is applied in the generated BOSS code.")

sel_gamKLamXiBar = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=3"
    nChrn     ">=3"
    nNet      "==0"
  }
  .select_photon {                      # Photon selection (radiative mode)
    tdc_emc_start     0                 # TDC window 0-14
    tdc_emc_end       14
    angle_to_track    10.0              # > 10 degrees to the nearest charged track
    energyThreshold_b 0.025             # 25 MeV in the barrel
    energyThreshold_e 0.050             # 50 MeV in the endcap
    nGam              ">=1"             # at least one photon
  }
  .pid(method: :probability) {
    prob_cut 0.001                      # CL > 0.001
    identify :proton, against: [:kaon, :pion]
    identify :kaon,   against: [:pion, :proton]
    nprp ">=1"
    nprm ">=1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:Lambda_bar, :pip]) {
    build_virtual_particle(:Xi_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # 4C kinematic fit: chi2 < 100, Lambda / anti-Lambda and anti-Xi mass windows.
  # If several photons are present the combination with the smallest chi2 is selected automatically.
  .kinematic_fit([:gamma, :km, :Lambda, :Lambda_bar, :Xi_bar]) {
    nominal
    constrain_four_momentum
    chi2_cut 100
    invariant_mass_of(:Lambda).within(1.110, 1.121)       # Lambda mass window (GeV/c^2)
    invariant_mass_of(:Lambda_bar).within(1.110, 1.121)   # anti-Lambda mass window (GeV/c^2)
    invariant_mass_of(:Xi_bar).within(1.315, 1.330)       # anti-Xi+ mass window (GeV/c^2)
  }

alg_gamKLamXiBar.with_decay_card(decay_card_gamKLamXiBar).apply(sel_gamKLamXiBar)

### Execute on datasets ###
root_files_KLamXiBar    = alg_KLamXiBar.execute_on([psip_data, psip_incMC, exMC_KLamXiBar])
root_files_gamKLamXiBar = alg_gamKLamXiBar.execute_on([psip_data, psip_incMC, exMC_gamKLamXiBar])