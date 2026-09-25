# BESIII DSL: 2405.09066v1 — First search for D*+ -> e+ nu_e and D*+ -> mu+ nu_mu
# e+e- collisions at sqrt(s) = 4.178–4.226 GeV, L = 6.32 fb^-1, BOSS 703
# ST+missing tag analysis: tag D*- via D0 hadronic modes, signal lepton + neutrino

# Datasets (BOSS 703, multi-energy D*+D*- production)
data_4180 = DatasetManager.real_data.find("703_4180")
data_4190 = DatasetManager.real_data.find("703_4190")
data_4200 = DatasetManager.real_data.find("703_4200")
data_4210 = DatasetManager.real_data.find("703_4210")
data_4220 = DatasetManager.real_data.find("703_4220")
data_4230 = DatasetManager.real_data.find("703_4230")

incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")

all_data = [data_4180, data_4190, data_4200, data_4210, data_4220, data_4230]
all_incMC = [incMC_4180, incMC_4190, incMC_4200, incMC_4210, incMC_4220, incMC_4230]

# Decay card — signal D*+ -> e+ nu_e (exclusive MC)
decay_card_enu = <<~DECAYCARD
    Decay psi(4260)
    1.000 D*+ D*- PHSP;
    Enddecay

    Decay D*+
    1.000 e+ nu_e PHSP;
    Enddecay

    Decay D*-
    1.000 anti-D0 pi- VSS;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay
    End
DECAYCARD

# Decay card — signal D*+ -> mu+ nu_mu (exclusive MC)
decay_card_munu = <<~DECAYCARD
    Decay psi(4260)
    1.000 D*+ D*- PHSP;
    Enddecay

    Decay D*+
    1.000 mu+ nu_mu PHSP;
    Enddecay

    Decay D*-
    1.000 anti-D0 pi- VSS;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay
    End
DECAYCARD

# Exclusive MC for signal processes
exMC_enu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "dstar_enu_exclusive_mc"
  config.related_dataset = data_4180
  config.events = 100000
  config.decay_card = decay_card_enu
  config.cross_section = :default
end

exMC_munu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "dstar_munu_exclusive_mc"
  config.related_dataset = data_4180
  config.events = 100000
  config.decay_card = decay_card_munu
  config.cross_section = :default
end

# =============================================
# Algorithm 1: D*+ -> e+ nu_e  (ST + missing nu_e)
# =============================================
alg_enu = TagAnalysis.new("DstarENu")
alg_enu.set_header(["DstarENuAlg/DstarENu.h"])
        .set_constant({ "ECMS" => [:double, 4.178] })
        .note(:tag_limitation, "Only D0 hadronic tag modes included (D0toKPi, D0toKPiPi0, D0toKPiPiPi); D+ tag modes (DptoKPiPi, DptoKKPi, DptoKPiPiPi0, DptoKsPi, DptoKsPiPi0, DptoKsPiPiPi) which also contribute to D*- ST via D- pi0 are not covered by this single-species TagAnalysis.")
        .note(:extra_photon_veto, "Maximum energy of extra photons not used by tag required < 0.3 GeV (no dedicated DSL method in TagAnalysis signal_side).")
        .note(:electron_pid, "Positron PID: L'(e)>0.001 and L'(e)/(L'(e)+L'(pi)+L'(K))>0.8, E/p>0.8. TagAnalysis uses fixed v1 lepton PID thresholds.")

alg_enu.tag_side(:D0) do |t|
  # D* -> D0bar pi-, D0bar -> hadrons; charm -1 tags the anti-D0
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

alg_enu.signal_side do |s|
  s.charged(ep: 1)        # one positron from D*+ -> e+ nu_e
  s.require_charge 1       # e+ charge = +1
  s.missing :nu_e           # massless neutrino (2-arg p4 AddMissTrack overload)
end

alg_enu.fit do |f|
  f.constrain_four_momentum   # 4C: tag D*- + e+ + nu_e = ecms_lab
  f.chi2_cut 200
end

alg_enu.with_decay_card(decay_card_enu).apply
alg_enu.execute_on(all_data + all_incMC + [exMC_enu])

# =============================================
# Algorithm 2: D*+ -> mu+ nu_mu  (ST + missing nu_mu)
# =============================================
alg_munu = TagAnalysis.new("DstarMuNu")
alg_munu.set_header(["DstarMuNuAlg/DstarMuNu.h"])
         .set_constant({ "ECMS" => [:double, 4.178] })
         .note(:tag_limitation, "Only D0 hadronic tag modes included; D+ tag modes (6 modes via D- pi0) are not covered.")
         .note(:extra_photon_veto, "Maximum energy of extra photons not used by tag required < 0.3 GeV.")
         .note(:muon_pid, "Muon PID: E_EMC in (0.0, 0.3) GeV, MUC hit depth requirements per cos(theta) and momentum bins. TagAnalysis uses fixed v1 muon PID thresholds.")

alg_munu.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

alg_munu.signal_side do |s|
  s.charged(mup: 1)        # one muon from D*+ -> mu+ nu_mu
  s.require_charge 1
  s.missing :nu_mu          # massless neutrino
end

alg_munu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_munu.with_decay_card(decay_card_munu).apply
alg_munu.execute_on(all_data + all_incMC + [exMC_munu])