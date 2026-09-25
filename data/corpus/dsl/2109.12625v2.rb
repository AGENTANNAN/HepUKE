# BESIII DSL: Search for di-muon decays of CP-odd light Higgs boson A0 in J/psi radiative decays
# arXiv: 2109.12625v2
# J/psi -> gamma A0, A0 -> mu+ mu-
# Uses 9.0 billion J/psi events at 3.097 GeV

### Dataset preparation ###
jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

decay_card = <<~DECAYCARD
    Decay J/psi
    1.000 gamma A0               SVP_HELAMP 1.0 0.0 1.0 0.0;
    Enddecay

    Decay A0
    1.000 mu+ mu-                PHSP;
    Enddecay

    End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_gamma_a0_mumu_exclusive_mc"
  config.related_dataset = jpsi_data
  config.events = 120_000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg = Algorithm.new("JpsiGammaA0MuMu")
alg.set_header(["JpsiGammaA0MuMuAlg/JpsiGammaA0MuMu.h"])
alg.set_constant({"ECMS" => [:double, 3.097]})

event_selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=1"
    nChrn ">=1"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=1"
  end
  # Lepton identification: muons identified via EMC energy, TOF timing, MUC depth
  .pid(method: :probability) do
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.5,
                                   treat_as_electron_if_energy_above: 0.6
    nlp ">=1"
    nlm ">=1"
  end
  # Vertex fit for mu+ mu- pair (A0 candidate)
  .kinematic_fit([:gamma, :lp, :lm]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg.note(:muon_pid, "Muon selection: E_cal^mu/p < 0.9c, 0.1 < E_cal^mu < 0.3 GeV, |Delta t_TOF| < 0.26 ns; MUC penetration depth > (-40+70*p) cm for p<1.1 GeV/c and >40 cm for p>1.1 GeV/c")
  .note(:vertex_fit, "mu+ mu- tracks required to originate from common vertex via vertex fit")
  .note(:kinematic_fit_4c, "4C kinematic fit with chi2_4C < 40; di-muon mass < 3.04 GeV/c^2 from fit")
  .note(:helicity_cut, "|cos(theta_hel_mu)| < 0.92 to suppress e+e- -> gamma mu+ mu- and J/psi -> mu+ mu- (gamma) backgrounds")
  .note(:reduced_mass, "m_red = sqrt(m(mu+mu-)^2 - 4 m_mu^2) used as fit variable; signal PDF modeled by sum of two Crystal Ball functions")
  .note(:background_modeling, "Non-peaking background (e+e- -> gamma mu+ mu- + J/psi -> mu+ mu- gamma) modeled by tanh(polynomial) at threshold, Chebyshev polynomial elsewhere; peaking background from rho/omega -> pi+pi- (GS function), gamma f with f=f2(1270)/f0(1500)/f0(1710) (double CB)")
  .note(:scan_search, "Search performed in 1 MeV/c^2 steps for m_A0 in [0.22, 1.5] GeV/c^2 and 2 MeV/c^2 steps in [1.5, 3.0] GeV/c^2; 2035 mass points total")
  .note(:upper_limits, "90% CL upper limits on B(J/psi->gamma A0)xB(A0->mu+mu-) in range (1.2-778.0)e-9 for 0.212 <= m_A0 <= 3.0 GeV/c^2; 6-7x improvement over previous BESIII result; slightly better than BaBar in low-mass region for tan(beta)=1")
  .with_decay_card(decay_card)
  .apply(event_selection)

root_files = alg.execute_on([jpsi_data, jpsi_incMC, exMC])