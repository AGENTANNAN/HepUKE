### Dataset preparation ###
# D+ -> omega mu+ nu_mu semileptonic decay at sqrt(s)=3.773 GeV, 2.93 fb-1
# Single-tag (ST) + missing neutrino using BOSS 712

data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Decay card: psi(3770) -> D+ D-, D+ -> omega mu+ nu_mu
# Generator-level card for signal MC
decay_card = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D-               VSS;
  Enddecay
  Decay D+
  1.0000 omega mu+ nu_mu     ISGW;
  Enddecay
  Decay D-
  1.0000 K+ pi- pi-          PHSP;
  Enddecay
  Decay omega
  1.0000 pi+ pi- pi0          PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma           PHSP;
  Enddecay
  End
DECAYCARD

# Exclusive MC for signal
exMC_sig = DatasetManager.create_exclusive_mc do |config|
  config.sample_name   = "sig_Dp_omega_munu"
  config.events        = 200_000
  config.decay_card    = decay_card
  config.related_dataset = data_3773
  config.cross_section = :default
end

### TagAnalysis: D+ -> omega mu+ nu_mu (ST + missing neutrino) ###
alg_omega_munu = TagAnalysis.new("Dp_omega_munu")
alg_omega_munu.set_header(["Dp_omega_munuAlg/Dp_omega_munu.h"])
               .set_constant({"ECMS" => [:double, 3.773]})
               .with_decay_card(decay_card)

# Tag side: single-tag D- (charm -1) with 6 hadronic modes
alg_omega_munu.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKSPi, :DptoKPiPiPi0,
          :DptoKSPiPi0, :DptoKSPiPiPi, :DptoKKPi
  t.charm -1
  t.window :deltaE, abs: 0.055      # broader window for pi0-containing tags
end

# Signal side: D+ -> omega mu+ nu_mu, omega -> pi+ pi- pi0
# pi+pi- from omega, mu+, pi0->gamma gamma, missing neutrino
alg_omega_munu.signal_side do |s|
  s.photons 2                        # pi0 -> gamma gamma
  s.charged(pip: 1, pim: 1, mup: 1)
  s.require_charge 1                 # D+ charge = +1
  s.missing :nu_mu                   # massless neutrino
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

# 4C kinematic fit: tag + signal + nu_mu = ecms_lab
alg_omega_munu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_omega_munu
  .note(:event_selection, "ST selected by M_BC in [1.863, 1.877] GeV and deltaE windows (broader for pi0 tags). Omega candidate: |M(pi+pi-pi0)-m_omega| < 0.025 GeV. KS0 veto: |M(pi+pi-)-m_KS0| > 0.015 GeV. Extra photon energy E_max < 0.15 GeV. M(omega mu+) < 1.5 GeV. Recoil mass veto: M_recoil(D-pi+pi+pi-) outside [0.45, 0.55] GeV.")
  .note(:dt_yield, "DT yield extracted from U_miss = E_miss - |p_miss| distribution via unbinned maximum likelihood fit. N_DT = 194 +/- 20. BF = N_DT / (N_ST_tot * epsilon_SL * B_omega * B_pi0).")
  .note(:form_factor, "Signal MC generated with ISGW model; form factor parameters r_V = 1.24, r_2 = 1.06. Systematic uncertainty from MC model evaluated by comparison with alternative ISGW2 model.")
  .note(:muon_pid, "Muon PID: CL_mu > 0.001, CL_mu > CL_e, CL_mu > CL_K, EMC deposited energy in (0.15, 0.25) GeV. Thresholds are BOSS-side defaults, not DSL-tunable.")
  .apply

# Execute on psi(3770) data
alg_omega_munu.execute_on([data_3773, incMC_3773, exMC_sig])