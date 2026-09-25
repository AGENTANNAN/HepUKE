# Paper 2312.08414v2: Born cross sections for e+e-→Λc+Λc(2595)-+c.c. and Λc+Λc(2625)-+c.c.
# Ordinary analysis at √s=4918.0 and 4950.9 MeV, BOSS 707

data_4918 = DatasetManager.load_real_data.find("707_4914")
data_4951 = DatasetManager.load_inclusive_mc.find("707_4946")
inc_mc_4918 = DatasetManager.load_inclusive_mc.find("707_4914")
inc_mc_4951 = DatasetManager.load_inclusive_mc.find("707_4946")

# Exclusive MC: e+e- → Λc+ Λc(2595)- and Λc+ Λc(2625)-
# Λc+ → p K- π+, excited Λc → Λc+ π π (generic)
excl_mc_4918 = DatasetManager.create_exclusive_mc do |c|
  c.decay_card = <<~DECAY
    Decay psi(4260)
    1.0 Lambda_c+ anti-Lambda_c(2595)- VSS;
    Decay Lambda_c+
    1.0 p K- pi+ PHSP;
    Decay anti-Lambda_c(2595)-
    1.0 anti-Lambda_c- pi+ pi- PHSP;
    Decay anti-Lambda_c-
    1.0 anti-p K+ pi- PHSP;
    End
  DECAY
  c.related_dataset data_4918
end

excl_mc_4951 = DatasetManager.create_exclusive_mc do |c|
  c.decay_card = <<~DECAY
    Decay psi(4260)
    1.0 Lambda_c+ anti-Lambda_c(2595)- VSS;
    Decay Lambda_c+
    1.0 p K- pi+ PHSP;
    Decay anti-Lambda_c(2595)-
    1.0 anti-Lambda_c- pi+ pi- PHSP;
    Decay anti-Lambda_c-
    1.0 anti-p K+ pi- PHSP;
    End
  DECAY
  c.related_dataset data_4951
end

# === Common selection for both energy points ===
# The selection is identical; the cross section and form-factor extraction is done in ROOT

# Algorithm for 4918.0 MeV
alg_4918 = Algorithm.new("LcLcStarXSection4918", version: "00-00-01")
alg_4918.set_header(["LcLcStarAlg/LcLcStar.h"])
  .set_constant({ "ECMS" => [:double, 4.918] })

sel_4918 = Selection.new("LcLcStarSel4918")

# Charged tracks: p, K-, π+ (minimum 3 tracks)
sel_4918.select_track do |t|
  t.nTot 3
  t.nChrp 0  # Λc+ = +1, but the recoil system is studied; require >=3 tracks
end

# PID: identify proton, kaon, pion by highest probability
sel_4918.pid(method: :probability) do |pid|
  pid.identify(:prp, :km, :pip, against: [:kp, :pim, :ep, :em, :mup, :mum])
  pid.prob_cut 0.001
end

# Build Λc+ candidate: p K- π+
sel_4918.build_virtual_particle(:Lambda_c, from: [:prp, :km, :pip])

# No kinematic fit needed — the recoiling mass M_rec+ is computed from the Λc+ candidate
# and used to identify Λc(2595)- and Λc(2625)- signals in ROOT

sel_4918.note(:Lambda_c_signal_region, "M(pK-π+) ∈ (2.27, 2.30) GeV/c2; M_rec+ > 2.55 GeV/c2")
sel_4918.note(:recoil_mass_analysis, "Unbinned ML fit to M_rec+ spectrum for signal yields of Λc(2595)- and Λc(2625)-")
sel_4918.note(:angular_distribution, "For Λc(2625): cosθ distribution with 6 bins at 4918.0 MeV, 8 bins at 4950.9 MeV")
sel_4918.note(:iterative_cross_section, "ISR correction 1+δ and efficiencies obtained iteratively using pQCD-motivated power-law line shape")
sel_4918.note(:background_components, "Backgrounds: e+e-→qq(continuum), e+e-→ΣcΣc, e+e-→Λc+Σcπ+c.c. — shapes from MC + ARGUS")
sel_4918.note(:track_quality, "Tracks: |cosθ|<0.93, Vz<10cm, Vxy<1cm (standard BESIII track selection)")
sel_4918.note(:excited_Lc_decay, "Λc(2595)+/Λc(2625)+ → Λc+π+π- and Λc+π0π0 with ratio 2:1 (isospin); ratio varied 1.5:1–5:1 for systematics")

alg_4918.apply(sel_4918)
alg_4918.execute_on([data_4918])

# Algorithm for 4950.9 MeV — same selection, different ECMS
alg_4951 = Algorithm.new("LcLcStarXSection4951", version: "00-00-01")
alg_4951.set_header(["LcLcStarAlg/LcLcStar.h"])
  .set_constant({ "ECMS" => [:double, 4.9509] })

sel_4951 = Selection.new("LcLcStarSel4951")

sel_4951.select_track do |t|
  t.nTot 3
  t.nChrp 0
end

sel_4951.pid(method: :probability) do |pid|
  pid.identify(:prp, :km, :pip, against: [:kp, :pim, :ep, :em, :mup, :mum])
  pid.prob_cut 0.001
end

sel_4951.build_virtual_particle(:Lambda_c, from: [:prp, :km, :pip])

sel_4951.note(:Lambda_c_signal_region, "M(pK-π+) ∈ (2.27, 2.30) GeV/c2; M_rec+ > 2.55 GeV/c2")
sel_4951.note(:recoil_mass_analysis, "Unbinned ML fit to M_rec+ spectrum for signal yields of Λc(2595)- and Λc(2625)-")
sel_4951.note(:angular_distribution, "For Λc(2625): cosθ distribution with 6 bins at 4918.0 MeV, 8 bins at 4950.9 MeV")
sel_4951.note(:iterative_cross_section, "ISR correction 1+δ and efficiencies obtained iteratively using pQCD-motivated power-law line shape")
sel_4951.note(:background_components, "Backgrounds: e+e-→qq(continuum), e+e-→ΣcΣc, e+e-→Λc+Σcπ+c.c. — shapes from MC + ARGUS")
sel_4951.note(:track_quality, "Tracks: |cosθ|<0.93, Vz<10cm, Vxy<1cm (standard BESIII track selection)")
sel_4951.note(:excited_Lc_decay, "Λc(2595)+/Λc(2625)+ → Λc+π+π- and Λc+π0π0 with ratio 2:1 (isospin); ratio varied 1.5:1–5:1 for systematics")

alg_4951.apply(sel_4951)
alg_4951.execute_on([data_4951])