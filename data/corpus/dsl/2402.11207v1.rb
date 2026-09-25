# Search for e+e- -> p p pi- dbar + c.c. at sqrt(s)=4.13-4.70 GeV
# BESIII: arXiv:2402.11207v1
# Partial reconstruction: deuteron (anti-deuteron) can be missed
# Four mutually exclusive event classes (3/4 tracks, charge-conjugate modes)

### Dataset — XYZ energy scan ###
data_points = [
  DatasetManager.real_data.find("705_4130"),
  DatasetManager.real_data.find("705_4160"),
  DatasetManager.real_data.find("703_4180"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4230"),
  DatasetManager.real_data.find("703_4260"),
  DatasetManager.real_data.find("705_4290"),
  DatasetManager.real_data.find("705_4315"),
  DatasetManager.real_data.find("705_4340"),
  DatasetManager.real_data.find("703_4360"),
  DatasetManager.real_data.find("705_4380"),
  DatasetManager.real_data.find("705_4400"),
  DatasetManager.real_data.find("703_4420"),
  DatasetManager.real_data.find("705_4440"),
  DatasetManager.real_data.find("706_4620"),
  DatasetManager.real_data.find("706_4640"),
  DatasetManager.real_data.find("706_4660"),
  DatasetManager.real_data.find("706_4680"),
  DatasetManager.real_data.find("706_4700"),
]
incMC_points = data_points.map { |d|
  begin
    DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}")
  rescue
    nil
  end
}.compact

# Decay card: e+e- -> p p pi- dbar (charge-conjugate included by c.c.)
decay_card_pppid = <<~DECAYCARD
  Decay psi(4260)
  1.000 p+ p+ pi- anti-d PHSP;
  Enddecay
  End
DECAYCARD

# Also need charge-conjugate: e+e- -> pbar pbar pi+ d
decay_card_ppbarpid = <<~DECAYCARD
  Decay psi(4260)
  1.000 anti-p- anti-p- pi+ d PHSP;
  Enddecay
  End
DECAYCARD

exMC_pppid = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_pppid"
  config.events        = 100_000
  config.decay_card    = decay_card_pppid
  config.cross_section = :default
end

exMC_ppbarpid = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_ppbarpid"
  config.events        = 100_000
  config.decay_card    = decay_card_ppbarpid
  config.cross_section = :default
end

### Algorithm: p p pi- + missing dbar (case 1: 3-track, deuteron missed) ###
alg_pppid = Algorithm.new("PpPidSearch")
alg_pppid.set_header(["PpPidSearchAlg/PpPidSearch.h"])

sel_pppid = Selection.new
sel_pppid.select_track do
  cos_theta 0.93
  Vz        10.0
  Vr        1.0
  nChrp     "==2"
  nChrn     "==1"
  nNet      "==1"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  identify :pion, against: [:kaon, :proton]
  nprp "==2"
  npim "==1"
end
.remove([:prp <= :chrgp, :pim <= :chrgn])
# Vertex fit: 2 protons + pi- from common vertex
.secondary_vertex_fit([:prp, :prp, :pim]) do
  build_virtual_particle(:vtx_pppi).by_minimizing_verfit_chi2
  remove_used_particle_from_candidate_list
end
# Partial reconstruction: reconstruct p p pi-, miss anti-d
.partial_miss([2]) do
  require_recoil_mass 1.80, 2.30
end

alg_pppid
  .note(:partial_reconstruction, "partial-reconstruction technique: deuteron can be missed; 4 mutually exclusive cases: case 1 (3trk, dbar missed), case 2 (4trk, dbar reconstructed), case 3 (3trk, d missed, c.c.), case 4 (4trk, d reconstructed, c.c.)")
  .note(:vertex_fit_chi2, "vertex fit chi2_VF < 70 for p p pi- vertex; vertex position constrained to (Rvx-0.12)^2 + (Rvy+0.14)^2 <= (0.5 cm)^2")
  .note(:pt_and_costheta_cut, "cases 1/3: pt(recoil) < 0.35 GeV/c when |cos(theta_r)| <= 0.93 to suppress background")
  .note(:md2_cut, "cases 2/4: deuteron mass squared from TOF: 3.2 < m_d^2 < 4.1 (GeV/c2)^2; cos(theta') > 0.9 between deuteron track and recoil direction")
  .note(:proton_veto, "cases 2/4: anti-proton veto in p p pi- dbar mode; proton veto in pbar pbar pi+ d mode to suppress p p pi- pbar nbar background")
  .note(:signal_region, "RM(p p pi-) and Rvz signal regions defined; no significant signal observed, upper limits set")
  .note(:geant4_deuteron, "deuteron not simulated in BESIII Geant4; efficiencies for cases 1/2 estimated from control samples; cross-efficiency assumptions validated with p pbar pi+ pi- control sample")
  .note(:no_significant_signal, "no significant signal observed; cross section upper limits set at 90% C.L. for each energy point")
  .with_decay_card(decay_card_pppid)
  .apply(sel_pppid)

### Charge-conjugate algorithm: pbar pbar pi+ + missing d (case 3) ###
alg_ppbarpid = Algorithm.new("PbarPbarPidSearch")
alg_ppbarpid.set_header(["PbarPbarPidSearchAlg/PbarPbarPidSearch.h"])

sel_ppbarpid = Selection.new
sel_ppbarpid.select_track do
  cos_theta 0.93
  Vz        10.0
  Vr        1.0
  nChrp     "==1"
  nChrn     "==2"
  nNet      "==-1"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  identify :pion, against: [:kaon, :proton]
  nprm "==2"
  npip "==1"
end
.remove([:prm <= :chrgn, :pip <= :chrgp])
.secondary_vertex_fit([:prm, :prm, :pip]) do
  build_virtual_particle(:vtx_ppbarpi).by_minimizing_verfit_chi2
  remove_used_particle_from_candidate_list
end
.partial_miss([2]) do
  require_recoil_mass 1.80, 2.30
end

alg_ppbarpid
  .note(:partial_reconstruction, "charge-conjugate channel: pbar pbar pi+ + missing d; case 3 of 4")
  .note(:vertex_fit, "vertex fit with chi2_VF < 70 for pbar pbar pi+ vertex")
  .with_decay_card(decay_card_ppbarpid)
  .apply(sel_ppbarpid)

alg_pppid.execute_on(data_points + incMC_points + exMC_pppid)
alg_ppbarpid.execute_on(data_points + incMC_points + exMC_ppbarpid)