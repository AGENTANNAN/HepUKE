# 1907.11370v3: Semileptonic D+ → K1(1270)^0 e+ ν_e at ψ(3770)

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

decay_card_Dp_K1enu = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D- PHSP;
  Enddecay
  Decay D+
  1.0000 K- pi+ pi0 e+ nu_e PHSP;
  Enddecay
  Decay D-
  1.0000 K+ pi- pi- PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

exMC_Dp_K1enu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_Dp_K1enu"
  config.related_dataset = psi3770_data
  config.events = 200_000
  config.decay_card = decay_card_Dp_K1enu
  config.cross_section = :default
end

# ====== TagAnalysis: ST D- + signal D+ → K-π+π0 e+ ν_e ======
alg = TagAnalysis.new("DpTagK1enu")
alg.set_header(["DpTagK1enuAlg/DpTagK1enu.h"])
    .set_constant({ "ECMS" => [:double, 3.773] })
    .with_decay_card(decay_card_Dp_K1enu)

# Tag side: D- via 6 hadronic modes
alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

# Signal side: K- π+ π0(→γγ) e+ + missing ν_e
alg.signal_side do |s|
  s.photons 2
  s.min_photon_angle 10.0
  s.charged(km: 1, pip: 1, ep: 1)
  s.missing :nu_e
end

# 5C fit: tag + K- + π+ + e+ + π0(mass-constrained) + ν_e = ecms_lab
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg.apply
alg.note(:lepton_pid, "Positron PID: CLe > 0.001 and CLe/(CLe+CLπ+CLK) > 0.8; E/p > 0.8 to suppress hadrons/muons")
    .note(:pi0_selection, "π0 momentum > 0.15 GeV/c; |cosθ_decay,π0| < 0.8; M(γγ) ∈ [0.115, 0.150] GeV; 1C mass-constrained fit applied")
    .note(:hadronic_veto, "M(K-π+π0 e+) < 1.78 GeV to suppress D+ → K-π+π+π0 hadronic bg")
    .note(:Kstar_veto, "U'_miss outside (-0.09, 0.03) GeV to veto K*(892)^0 e+ν with fake π0")
    .note(:fsr_recovery, "FSR/bremsstrahlung photons within 5° of positron added to e+ 4-momentum")
    .note(:signal_model, "D+→K1eν simulated with ISGW2 model; K1 resonance shape via relativistic Breit-Wigner")
    .note(:twoD_fit, "2D unbinned ML fit to M(K-π+π0) vs U_miss; signal from MC shape; bg from inclusive MC shape")

alg.execute_on([psi3770_data, psi3770_incMC, exMC_Dp_K1enu])