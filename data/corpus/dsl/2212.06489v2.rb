# 2212.06489v2: CP-even fraction F+ in D0 -> K+K-pi+pi- at psi(3770)
# TagAnalysis: D-tag with CP-eigenstate tags, signal D0 -> K+K-pi+pi-
# Dataset: 712_3773 (psi(3770))

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card: psi(3770) -> D0 anti-D0
# D0 -> K+ K- pi+ pi-  (signal)
# anti-D0 -> [CP-tag decay modes] (tag side)
# Conventions: D0 charge conjugate implied throughout
decay_card = <<~DECAY
Decay psi(3770)
  1.000 D0 anti-D0 VSS;
Enddecay
Decay D0
  1.000 K+ K- pi+ pi- PHSP;
Enddecay
Decay anti-D0
  1.000 K+ pi- PHSP;
Enddecay
End
DECAY

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_D0_KKpipi"
  config.related_dataset = psi3770_data
  config.events          = 1_000_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

alg = TagAnalysis.new("D0TagCPFrac", "00-00-01")
alg.set_header(["D0TagCPFrac/D0TagCPFrac.h"])
alg.set_constant({ "ECMS" => [:double, 3.773] })
alg.with_decay_card(decay_card)

# Tag side: anti-D0 tagged in CP eigenstates
# Paper Table I tag modes:
#   CP even: K+K-, pi+pi-, KS0pi0pi0, pi+pi-pi0, KL0pi0
#   CP odd: KS0pi0, KS0 eta(gg), KS0 eta'(pipi eta), KS0 eta'(rho0 gamma), KS0 omega
#   Mixed CP: KS0pi+pi-, KL0pi+pi-
alg.tag_side(:D0) do |t|
  # Standard DTagAlg D0 modes available in the DSL:
  t.modes :D0toKPiPiPi       # K pi pi pi (closest to K+K-pi+pi- but different)
  t.modes :D0toKsPiPi         # KS0 pi+ pi-  (Mixed CP tag)
  t.charm -1                  # tag anti-D0
end

alg.note(:tag_mode_unavailable,
  "Paper tag modes listed in Table I are NOT in the frozen BOSS 7.0.6/7.1.2 EvtRecDTag enum. " \
  "CP-even: K+K-, pi+pi-, KS0pi0pi0, pi+pi-pi0, KL0pi0. " \
  "CP-odd: KS0pi0, KS0 eta(gammagamma), KS0 eta'(pipi eta), KS0 eta'(rho0 gamma), KS0 omega. " \
  "Mixed CP: KS0pi+pi- (available as :D0toKsPiPi), KL0pi+pi-. " \
  "Only two standard D0 modes declared above; full CP-tag mode list requires custom EvtRecDTag.h fork. " \
  "Paper also has partial reconstruction for KL0 modes (missing-mass-squared fits) and Kalman kinematic fit for KS,L0pi+pi- tags."
)

# Signal side: D0 -> K+ K- pi+ pi-
# Four charged tracks from signal D
alg.signal_side do |s|
  s.charged(km: 1, kp: 1, pip: 1, pim: 1)
  s.require_charge 0            # K+K-pi+pi- is neutral total
end

alg.note(:signal_side_selection,
  "Signal side selection: pi+pi- pair required to originate from vertex within 2x vertex resolution " \
  "from IP to suppress D->KS0KK background. pi+pi- invariant mass outside [477, 507] MeV/c^2 " \
  "(KS0 veto). Charged kaons/pions identified by L(K)>L(pi) and L(pi)>L(K)."
)

alg.note(:tag_side_selection,
  "Tag side: KS0 reconstructed from pi+pi- (|Vz|<20cm, secondary vertex fit, mass within 12 MeV/c^2, " \
  "decay length > 2x vertex resolution). pi0/eta from di-photon: M(gammagamma) in [115,150] and " \
  "[480,580] MeV/c^2. eta' reconstructed via pi+pi-eta ([940,976] MeV/c^2) " \
  "and rho0 gamma ([940,970] MeV/c^2, M(pipi) in [626,924] MeV/c^2). " \
  "pi+pi-pi0: KS0 veto (M(pipi) > 18 MeV/c^2 from KS0 mass). " \
  "Fully reconstructed tags: DeltaE = E_D - sqrt(s)/2 within 3sigma. " \
  "KL0 modes: partial reconstruction using missing momentum. " \
  "Kalman kinematic fit for KS,L0pi+pi- tags constraining KS,L0 and D masses."
)

alg.note(:yield_extraction,
  "ST yields from ML fit to M_BC = sqrt(E_beam^2 - |sum p_i|^2). " \
  "Signal shape: MC shape convoluted with Gaussian. Bkg: ARGUS function. " \
  "DT yields from M_BC fits (fully reconstructed) or M_miss^2 fits (partially reconstructed). " \
  "KS0omega: sPlot technique on M_BC then pi+pi-pi0 mass fit. " \
  "Binned phase-space analysis for KS0pi+pi- and KL0pi+pi- using 8-bin scheme."
)

alg.note(:cp_fraction_fit,
  "F+ extracted from simultaneous ML fit to ST and DT yields of CP tags " \
  "using Eq.(3) and KS,L0pipi tags using Eq.(4). " \
  "External inputs: F+(pipipi0)=0.973, Ki/ci from CLEO+BESIII combined results."
)

# Fit: 4-momentum conservation
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg.apply
alg.execute_on([psi3770_data, psi3770_incMC, exMC])