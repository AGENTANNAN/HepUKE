# Paper 2004.01394v2: chi_cJ -> Sigma- anti-Sigma+ BF measurement at psi(3686)
# Decay: psi(3686) -> gamma chi_cJ, chi_cJ -> Sigma- anti-Sigma+, Sigma- -> n pi-, anti-Sigma+ -> anti-n pi+
# Data: psi(2S) at 3.686 GeV
# Special: anti-neutron selected via EMC showers, neutron treated as MISSING in kinematic fit

decay_card = <<~DECAYCARD
Decay psi(2S)
1.0 gamma chi_c0 PHSP;
Enddecay
Decay psi(2S)
1.0 gamma chi_c1 PHSP;
Enddecay
Decay psi(2S)
1.0 gamma chi_c2 PHSP;
Enddecay
Decay chi_c0
1.0 Sigma- anti-Sigma+ PHSP;
Enddecay
Decay chi_c1
1.0 Sigma- anti-Sigma+ PHSP;
Enddecay
Decay chi_c2
1.0 Sigma- anti-Sigma+ PHSP;
Enddecay
Decay Sigma-
1.0 n0 pi- PHSP;
Enddecay
Decay anti-Sigma+
1.0 anti-n0 pi+ PHSP;
Enddecay
End
DECAYCARD

psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "chi_cJ_Sigma_SigmaBar_exclusive_mc"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

algorithm = Algorithm.new("chi_cJ_to_Sigma_SigmaBar")
algorithm.set_header(["chi_cJ_to_Sigma_SigmaBarAlg/chi_cJ_to_Sigma_SigmaBar.h"])
algorithm.set_constant({ "ECMS" => [:double, 3.686] })

# Anti-neutron is selected from the most energetic EMC shower with dedicated criteria
# not expressible in DSL (energy 0.2-2.0 GeV, second moment S > 20 cm^2,
# >20 EMC hits in 40-degree cone). Neutron is treated as a MISSING particle.
# The kinematic fit: psi(3686) -> gamma + anti-n + pi+ + pi- with n MISSING,
# anti-n pi+ invariant mass constrained to anti-Sigma+ nominal mass.
# chi2 < 20; best combination by minimum chi2.
algorithm.note(:anti_neutron_selection,
  "Anti-neutron selected from most energetic EMC shower: deposited energy 0.2-2.0 GeV, " \
  "second moment S > 20 cm^2, >20 EMC hits in 40-degree cone. Purity > 98% from signal MC."
)

algorithm.note(:kinematic_fit,
  "Kinematic fit: psi(3686) -> gamma + anti-n + pi+ + pi- with n MISSING. " \
  "anti-n pi+ invariant mass constrained to anti-Sigma+ nominal mass. " \
  "chi2 < 20. Best combination by minimum chi2. " \
  "Cannot be expressed in DSL as anti-n and n are non-standard particle types " \
  "requiring custom EMC selection and missing-particle handling outside DSL framework."
)

algorithm.note(:background_vetoes,
  "Pre-kinematic-fit vetoes (applied in C++ code): " \
  "1. pi0 veto: reject event if any two gamma within |M_gg - m_pi0| < 12 MeV. " \
  "2. K_S0 veto: |M(pi+pi-) - m(K_S0)| > 10 MeV. " \
  "3. J/psi veto: |M_rec(pi+pi-) - m(J/psi)| > 10 MeV. " \
  "4. Competing-hypothesis veto: reject if chi2(Sigma- anti-Sigma+) < chi2(gamma Sigma- anti-Sigma+). " \
  "Radiative photon: opening angle > 40 deg from anti-n, > 10 deg from tracks, energy > 80 MeV."
)

# Minimal selection chain: track/photon selection + PID for pions
event_selection = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 30.0
    Vr 10.0
    nChrp "==1"
    nChrn "==1"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip "==1"
    npim "==1"
  }

# The kinematic fit cannot be expressed in DSL because anti-n and n are
# not standard particle types. The entire fit with miss_track_of(:n)
# and anti-n pi+ mass constraint must be implemented in generated C++ code.

algorithm.with_decay_card(decay_card).apply(event_selection)
algorithm.execute_on([psip_data, psip_incMC, exMC_signal])