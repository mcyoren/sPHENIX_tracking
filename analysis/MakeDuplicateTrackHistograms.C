// MakeDuplicateTrackHistograms.C
//
// Example:
// root -l -b -q 'MakeDuplicateTrackHistograms.C("/path/to/files","*.root","duplicate_tracks.root")'
//
// The macro compares every unique pair of selected tracks inside the same event.

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TMath.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace
{
  struct Track
  {
    double px = 0.;
    double py = 0.;
    double pz = 0.;
    double pt = 0.;
    double charge = 0.;
    double eta = 0.;
    double theta = 0.;
    double phi = 0.;
    double pca_z = 0.;
    double rDCA_zero = 0.;
    unsigned int ntpc_clusters = 0;
    int side = -1;
    unsigned int track_id = 0;
  };

  double deltaPhi(const double phi1, const double phi2)
  {
    return std::atan2(
      std::sin(phi1 - phi2),
      std::cos(phi1 - phi2)
    );
  }

  void compareTracks(
    const std::vector<Track>& tracks,
    TH3D* h_dpx_dpy_dpz,
    TH3D* h_dphi_dpt_dtheta,
    TH2D* h_dphi_dpt,
    TH2D* h_dphi_relDpt,
    TH2D* h_deltaP_relDpt,
    TH2D* h_deltaP_dtheta,
    TH1D* h_deltaP,
    TH1D* h_deltaPhi,
    TH1D* h_deltaPt,
    TH1D* h_relDeltaPt,
    TH1D* h_deltaTheta,
    TH1D* h_duplicateScore,
    TH1D* h_closePairsPerEvent,
    TH1D* h_trackPairsPerEvent)
  {
    long long nPairs = 0;
    long long nClosePairs = 0;

    for (std::size_t i = 0; i < tracks.size(); ++i)
    {
      for (std::size_t j = i + 1; j < tracks.size(); ++j)
      {
        const Track& a = tracks[i];
        const Track& b = tracks[j];

        ++nPairs;

        const double dpx = a.px - b.px;
        const double dpy = a.py - b.py;
        const double dpz = a.pz - b.pz;

        const double dpt = a.pt - b.pt;
        const double absDpt = std::abs(dpt);

        const double meanPt = 0.5 * (a.pt + b.pt);
        const double relDpt =
          meanPt > 0. ? absDpt / meanPt : 999.;

        const double dphi = deltaPhi(a.phi, b.phi);
        const double absDphi = std::abs(dphi);

        const double dtheta = a.theta - b.theta;
        const double absDtheta = std::abs(dtheta);

        const double deltaP =
          std::sqrt(dpx*dpx + dpy*dpy + dpz*dpz);

        // Dimensionless combined closeness measure.
        // Values near zero are the most duplicate-like.
        const double duplicateScore =
          std::sqrt(
            std::pow(absDphi / 0.002, 2) +
            std::pow(relDpt / 0.005, 2) +
            std::pow(absDtheta / 0.002, 2)
          );

        h_dpx_dpy_dpz->Fill(dpx, dpy, dpz);
        h_dphi_dpt_dtheta->Fill(dphi, dpt, dtheta);

        h_dphi_dpt->Fill(dphi, dpt);
        h_dphi_relDpt->Fill(dphi, relDpt);
        h_deltaP_relDpt->Fill(deltaP, relDpt);
        h_deltaP_dtheta->Fill(deltaP, dtheta);

        h_deltaP->Fill(deltaP);
        h_deltaPhi->Fill(dphi);
        h_deltaPt->Fill(dpt);
        h_relDeltaPt->Fill(relDpt);
        h_deltaTheta->Fill(dtheta);
        h_duplicateScore->Fill(duplicateScore);

        // Preliminary duplicate-candidate definition.
        // It is intentionally tight but can be changed after looking at the plots.
        const bool duplicateCandidate =
          absDphi < 0.002 &&
          relDpt < 0.005 &&
          absDtheta < 0.002;

        if (duplicateCandidate)
        {
          ++nClosePairs;
        }
      }
    }

    h_trackPairsPerEvent->Fill(nPairs);
    h_closePairsPerEvent->Fill(nClosePairs);
  }
}

void MakeDuplicateTrackHistograms(
  const char* inputDir = ".",
  const char* filePattern = "*.root",
  const char* outputName = "duplicate_track_histograms.root",
  const char* treeName = "residuals")
{
  const TString inputPattern =
    TString::Format("%s/%s", inputDir, filePattern);

  TChain chain(treeName);
  const int nFiles = chain.Add(inputPattern);

  if (nFiles <= 0)
  {
    std::cerr << "ERROR: no files matched "
              << inputPattern << std::endl;
    return;
  }

  if (!chain.GetBranch("event"))
  {
    std::cerr << "ERROR: tree does not contain an event branch."
              << std::endl;
    return;
  }

  // Change these types if ROOT reports a branch-type mismatch.
  UInt_t event = 0;
  UInt_t final_track_id = 0;
  UInt_t ntpc_clusters = 0;

  Int_t side = -1;

  Double_t px = 0.;
  Double_t py = 0.;
  Double_t pz = 0.;
  Double_t pt = 0.;
  Double_t eta = 0.;
  Double_t charge = 0.;
  Double_t pca_z = 0.;
  Double_t rDCA_zero = 0.;

  chain.SetBranchAddress("event", &event);
  chain.SetBranchAddress("px", &px);
  chain.SetBranchAddress("py", &py);
  chain.SetBranchAddress("pz", &pz);
  chain.SetBranchAddress("pt", &pt);
  chain.SetBranchAddress("eta", &eta);
  chain.SetBranchAddress("charge", &charge);
  chain.SetBranchAddress("side", &side);
  chain.SetBranchAddress("ntpc_clusters", &ntpc_clusters);
  chain.SetBranchAddress("pca_z", &pca_z);
  chain.SetBranchAddress("rDCA_zero", &rDCA_zero);

  const bool hasTrackId =
    chain.GetBranch("final_track_id") != nullptr;

  if (hasTrackId)
  {
    chain.SetBranchAddress(
      "final_track_id",
      &final_track_id
    );
  }

  TFile output(outputName, "RECREATE");

  TH3D* h_dpx_dpy_dpz =
    new TH3D(
      "h_dpx_dpy_dpz",
      "Same-event track pairs;"
      "#Deltap_{x} [GeV/c];"
      "#Deltap_{y} [GeV/c];"
      "#Deltap_{z} [GeV/c]",
      200, -0.20, 0.20,
      200, -0.20, 0.20,
      200, -0.20, 0.20
    );

  TH3D* h_dphi_dpt_dtheta =
    new TH3D(
      "h_dphi_dpt_dtheta",
      "Same-event track pairs;"
      "#Delta#phi [rad];"
      "#Deltap_{T} [GeV/c];"
      "#Delta#theta [rad]",
      200, -0.05, 0.05,
      200, -0.20, 0.20,
      200, -0.05, 0.05
    );

  TH2D* h_dphi_dpt =
    new TH2D(
      "h_dphi_dpt",
      "Same-event track pairs;"
      "#Delta#phi [rad];"
      "#Deltap_{T} [GeV/c]",
      400, -0.05, 0.05,
      400, -0.20, 0.20
    );

  TH2D* h_dphi_relDpt =
    new TH2D(
      "h_dphi_relDpt",
      "Same-event track pairs;"
      "#Delta#phi [rad];"
      "|#Deltap_{T}|/#LTp_{T}#GT",
      400, -0.05, 0.05,
      300, 0., 0.15
    );

  TH2D* h_deltaP_relDpt =
    new TH2D(
      "h_deltaP_relDpt",
      "Same-event track pairs;"
      "|#Delta#vec{p}| [GeV/c];"
      "|#Deltap_{T}|/#LTp_{T}#GT",
      300, 0., 0.30,
      300, 0., 0.15
    );

  TH2D* h_deltaP_dtheta =
    new TH2D(
      "h_deltaP_dtheta",
      "Same-event track pairs;"
      "|#Delta#vec{p}| [GeV/c];"
      "#Delta#theta [rad]",
      300, 0., 0.30,
      400, -0.05, 0.05
    );

  TH1D* h_deltaP =
    new TH1D(
      "h_deltaP",
      "Same-event track pairs;"
      "|#Delta#vec{p}| [GeV/c];pairs",
      500, 0., 0.50
    );

  TH1D* h_deltaPhi =
    new TH1D(
      "h_deltaPhi",
      "Same-event track pairs;"
      "#Delta#phi [rad];pairs",
      1000, -0.10, 0.10
    );

  TH1D* h_deltaPt =
    new TH1D(
      "h_deltaPt",
      "Same-event track pairs;"
      "#Deltap_{T} [GeV/c];pairs",
      1000, -0.50, 0.50
    );

  TH1D* h_relDeltaPt =
    new TH1D(
      "h_relDeltaPt",
      "Same-event track pairs;"
      "|#Deltap_{T}|/#LTp_{T}#GT;pairs",
      500, 0., 0.25
    );

  TH1D* h_deltaTheta =
    new TH1D(
      "h_deltaTheta",
      "Same-event track pairs;"
      "#Delta#theta [rad];pairs",
      1000, -0.10, 0.10
    );

  TH1D* h_duplicateScore =
    new TH1D(
      "h_duplicateScore",
      "Combined duplicate-track score;"
      "duplicate score;pairs",
      500, 0., 50.
    );

  TH1D* h_trackMultiplicity =
    new TH1D(
      "h_trackMultiplicity",
      "Selected tracks per event;"
      "tracks per event;events",
      200, -0.5, 199.5
    );

  TH1D* h_trackPairsPerEvent =
    new TH1D(
      "h_trackPairsPerEvent",
      "Track pairs per event;"
      "pairs per event;events",
      500, -0.5, 4999.5
    );

  TH1D* h_closePairsPerEvent =
    new TH1D(
      "h_closePairsPerEvent",
      "Duplicate-like pairs per event;"
      "duplicate-like pairs per event;events",
      50, -0.5, 49.5
    );

  std::vector<Track> eventTracks;

  ULong64_t currentEvent = 0;
  int currentTreeNumber = -1;
  bool firstAcceptedTrack = true;

  const Long64_t nEntries = chain.GetEntries();

  for (Long64_t entry = 0; entry < nEntries; ++entry)
  {
    chain.GetEntry(entry);

    const int treeNumber = chain.GetTreeNumber();

    // Important for a TChain:
    // event numbers may restart in every input file.
    const bool newEvent =
      !firstAcceptedTrack &&
      (treeNumber != currentTreeNumber ||
       event != currentEvent);

    if (newEvent)
    {
      h_trackMultiplicity->Fill(eventTracks.size());

      compareTracks(
        eventTracks,
        h_dpx_dpy_dpz,
        h_dphi_dpt_dtheta,
        h_dphi_dpt,
        h_dphi_relDpt,
        h_deltaP_relDpt,
        h_deltaP_dtheta,
        h_deltaP,
        h_deltaPhi,
        h_deltaPt,
        h_relDeltaPt,
        h_deltaTheta,
        h_duplicateScore,
        h_closePairsPerEvent,
        h_trackPairsPerEvent
      );

      eventTracks.clear();
    }

    if (firstAcceptedTrack || newEvent)
    {
      currentEvent = event;
      currentTreeNumber = treeNumber;
      firstAcceptedTrack = false;
    }

    if (!std::isfinite(px) ||
        !std::isfinite(py) ||
        !std::isfinite(pz) ||
        !std::isfinite(pt))
    {
      continue;
    }

    // Starting track selection.
    // Do not cut rDCA_zero here because duplicate tracks may share
    // the same bad DCA and are still worth detecting.
    if (ntpc_clusters < 20) continue;
    if (std::abs(pca_z) > 10.) continue;
    if (pt < 0.2) continue;

    Track track;

    track.px = px;
    track.py = py;
    track.pz = pz;
    track.pt = pt;
    track.eta = eta;
    track.charge = charge;
    track.side = side;
    track.pca_z = pca_z;
    track.rDCA_zero = rDCA_zero;
    track.ntpc_clusters = ntpc_clusters;
    track.track_id = hasTrackId ? final_track_id : 0;

    track.phi = std::atan2(py, px);

    const double momentum =
      std::sqrt(px*px + py*py + pz*pz);

    track.theta =
      momentum > 0.
      ? std::acos(
          std::clamp(pz / momentum, -1., 1.)
        )
      : 0.;

    eventTracks.push_back(track);
  }

  // Process the last event.
  if (!eventTracks.empty())
  {
    h_trackMultiplicity->Fill(eventTracks.size());

    compareTracks(
      eventTracks,
      h_dpx_dpy_dpz,
      h_dphi_dpt_dtheta,
      h_dphi_dpt,
      h_dphi_relDpt,
      h_deltaP_relDpt,
      h_deltaP_dtheta,
      h_deltaP,
      h_deltaPhi,
      h_deltaPt,
      h_relDeltaPt,
      h_deltaTheta,
      h_duplicateScore,
      h_closePairsPerEvent,
      h_trackPairsPerEvent
    );
  }

  output.Write();
  output.Close();

  std::cout
    << "Wrote duplicate-track histograms to "
    << outputName << std::endl;
}
