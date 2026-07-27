// Make2DClusterCorrectionVotes_v2.C
//
// First non-factorized test of the transverse 2D correction voting idea.
// Only clusters with |z_cluster| < maxAbsClusterZ are used.
// Histograms are organized as
//   side / charge / pT bin / cluster-phi bin
// and contain Delta r and r*Delta phi versus TPC layer.
//
// Important geometry:
// The closest point on a reference track gives only the component of the
// correction normal to that track. It does NOT determine the component along
// the local track direction. Therefore each cluster also fills a LINE of
// allowed (Delta r, r*Delta phi) corrections:
//
//   n_r * Delta r + n_phi * (r Delta phi) = d_normal.
//
// The intersection / maximum of many such lines is the actual 2D vote.
//
// Example:
// root -l -b -q 'Make2DClusterCorrectionVotes_v2.C("/path","*.root","votes_v2.root")'

#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THnSparse.h>
#include <TMath.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace
{
  constexpr double kTwoPi = 2.0 * TMath::Pi();

  struct XYPoint
  {
    double x = 0.;
    double y = 0.;
  };

  struct Circle
  {
    double xc = 0.;
    double yc = 0.;
    double radius = 0.;
    bool valid = false;
  };

  struct GroupHistograms
  {
    TH2D* h_deltaR_vs_layer = nullptr;
    TH2D* h_deltaRPhi_vs_layer = nullptr;
    TH2D* h_normalRadial_vs_layer = nullptr;
    TH2D* h_normalTangential_vs_layer = nullptr;
    TH2D* h_normalDistance_vs_layer = nullptr;
    TH2D* h_closestPoint_deltaR_deltaRPhi = nullptr;

    // Sparse axes: layer, Delta r, r Delta phi.
    // Every cluster fills a sampled line of allowed 2D corrections.
    THnSparseD* h_lineVotes_layer_deltaR_deltaRPhi = nullptr;
  };

  double wrapPhi(double phi)
  {
    while (phi >= TMath::Pi()) phi -= kTwoPi;
    while (phi < -TMath::Pi()) phi += kTwoPi;
    return phi;
  }

  double distance2(const XYPoint& a, const XYPoint& b)
  {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    return dx * dx + dy * dy;
  }

  Circle buildTrackCircle(
    const double pcaX,
    const double pcaY,
    const double px,
    const double py,
    const double pt,
    const double charge,
    const double bFieldTesla)
  {
    Circle circle;
    const double denom = 0.003 * charge * bFieldTesla;

    if (!std::isfinite(denom) || std::abs(denom) < 1e-12 ||
        !std::isfinite(pt) || pt <= 0.)
    {
      return circle;
    }

    circle.xc = pcaX + py / denom;
    circle.yc = pcaY - px / denom;
    circle.radius = std::abs(pt / denom);
    circle.valid = std::isfinite(circle.xc) &&
                   std::isfinite(circle.yc) &&
                   std::isfinite(circle.radius) &&
                   circle.radius > 0.;
    return circle;
  }

  std::vector<XYPoint> intersectCircleWithOriginRadius(
    const Circle& circle,
    const double detectorRadius)
  {
    std::vector<XYPoint> result;
    if (!circle.valid || detectorRadius <= 0.) return result;

    const double d = std::hypot(circle.xc, circle.yc);
    if (d < 1e-12) return result;

    const double r0 = detectorRadius;
    const double r1 = circle.radius;
    if (d > r0 + r1 || d < std::abs(r0 - r1)) return result;

    const double a = (r0 * r0 - r1 * r1 + d * d) / (2. * d);
    double h2 = r0 * r0 - a * a;
    if (h2 < -1e-9) return result;
    h2 = std::max(0.0, h2);

    const double h = std::sqrt(h2);
    const double x2 = a * circle.xc / d;
    const double y2 = a * circle.yc / d;
    const double rx = -circle.yc * h / d;
    const double ry =  circle.xc * h / d;

    result.push_back({x2 + rx, y2 + ry});
    if (h > 1e-10) result.push_back({x2 - rx, y2 - ry});
    return result;
  }

  double circleAnchorCost(
    const Circle& original,
    const Circle& candidate,
    const std::vector<double>& anchorRadii)
  {
    double cost = 0.;
    int used = 0;

    for (const double r : anchorRadii)
    {
      const auto oldPoints = intersectCircleWithOriginRadius(original, r);
      const auto newPoints = intersectCircleWithOriginRadius(candidate, r);
      if (oldPoints.empty() || newPoints.empty()) continue;

      double best = std::numeric_limits<double>::max();
      for (const auto& oldPoint : oldPoints)
      {
        for (const auto& newPoint : newPoints)
        {
          best = std::min(best, distance2(oldPoint, newPoint));
        }
      }

      if (std::isfinite(best))
      {
        cost += best;
        ++used;
      }
    }

    if (used < 2) return std::numeric_limits<double>::max();
    return cost / static_cast<double>(used);
  }

  Circle constrainCircleToBeamAndR2(
    const Circle& original,
    const double beamX,
    const double beamY,
    const std::vector<double>& anchorRadii,
    const int coarseScanSteps = 720,
    const int refinementSteps = 4)
  {
    Circle best;
    if (!original.valid) return best;

    double bestTheta = 0.;
    double bestCost = std::numeric_limits<double>::max();

    for (int i = 0; i < coarseScanSteps; ++i)
    {
      const double theta = kTwoPi * static_cast<double>(i) /
                           static_cast<double>(coarseScanSteps);
      Circle candidate;
      candidate.xc = beamX + original.radius * std::cos(theta);
      candidate.yc = beamY + original.radius * std::sin(theta);
      candidate.radius = original.radius;
      candidate.valid = true;

      const double cost = circleAnchorCost(original, candidate, anchorRadii);
      if (cost < bestCost)
      {
        bestCost = cost;
        bestTheta = theta;
        best = candidate;
      }
    }

    if (!best.valid) return best;

    double halfWidth = kTwoPi / static_cast<double>(coarseScanSteps);
    for (int level = 0; level < refinementSteps; ++level)
    {
      const int nLocal = 80;
      double localBestTheta = bestTheta;
      double localBestCost = bestCost;

      for (int i = 0; i <= nLocal; ++i)
      {
        const double frac = static_cast<double>(i) / nLocal;
        const double theta = bestTheta - halfWidth + 2. * halfWidth * frac;

        Circle candidate;
        candidate.xc = beamX + original.radius * std::cos(theta);
        candidate.yc = beamY + original.radius * std::sin(theta);
        candidate.radius = original.radius;
        candidate.valid = true;

        const double cost = circleAnchorCost(original, candidate, anchorRadii);
        if (cost < localBestCost)
        {
          localBestCost = cost;
          localBestTheta = theta;
          best = candidate;
        }
      }

      bestTheta = localBestTheta;
      bestCost = localBestCost;
      halfWidth *= 0.2;
    }

    return best;
  }

  XYPoint closestPointOnCircle(const XYPoint& point, const Circle& circle)
  {
    const double dx = point.x - circle.xc;
    const double dy = point.y - circle.yc;
    const double norm = std::hypot(dx, dy);

    if (!circle.valid || norm < 1e-12)
    {
      const double nan = std::numeric_limits<double>::quiet_NaN();
      return {nan, nan};
    }

    return {
      circle.xc + circle.radius * dx / norm,
      circle.yc + circle.radius * dy / norm
    };
  }

  int ptBin(const double pt)
  {
    if (pt >= 0.2 && pt < 0.7) return 0;
    if (pt >= 0.7 && pt < 1.5) return 1;
    if (pt >= 1.5 && pt < 3.0) return 2;
    if (pt >= 3.0) return 3;
    return -1;
  }

  std::string ptBinName(const int bin)
  {
    if (bin == 0) return "pt_0p2_0p7";
    if (bin == 1) return "pt_0p7_1p5";
    if (bin == 2) return "pt_1p5_3";
    if (bin == 3) return "pt_gt3";
    return "pt_invalid";
  }

  int phiBin(const double phi, const int nPhiBins)
  {
    const double wrapped = wrapPhi(phi);
    int bin = static_cast<int>((wrapped + TMath::Pi()) / kTwoPi * nPhiBins);
    if (bin < 0) bin = 0;
    if (bin >= nPhiBins) bin = nPhiBins - 1;
    return bin;
  }

  std::string groupName(
    const int side,
    const double charge,
    const int iPt,
    const int iPhi)
  {
    return TString::Format(
      "side%d/%s/%s/phi_%02d",
      side,
      charge > 0. ? "qplus" : "qminus",
      ptBinName(iPt).c_str(),
      iPhi).Data();
  }

  GroupHistograms bookGroup(
    TDirectory* dir,
    const std::string& label,
    const double maxAbsVoteCm)
  {
    dir->cd();
    GroupHistograms h;
    const TString suffix = TString::Format(" [%s, |z_{cluster}|<15 cm]", label.c_str());

    h.h_deltaR_vs_layer = new TH2D(
      "h_deltaR_vs_layer",
      "Closest-point #Delta r versus layer" + suffix +
      ";TPC layer;#Delta r [cm]",
      48, 6.5, 54.5, 160, -maxAbsVoteCm, maxAbsVoteCm);

    h.h_deltaRPhi_vs_layer = new TH2D(
      "h_deltaRPhi_vs_layer",
      "Closest-point r#Delta#phi versus layer" + suffix +
      ";TPC layer;r#Delta#phi [cm]",
      48, 6.5, 54.5, 160, -maxAbsVoteCm, maxAbsVoteCm);

    h.h_normalRadial_vs_layer = new TH2D(
      "h_normalRadial_vs_layer",
      "Radial component of track normal versus layer" + suffix +
      ";TPC layer;n_{r}",
      48, 6.5, 54.5, 120, -1.2, 1.2);

    h.h_normalTangential_vs_layer = new TH2D(
      "h_normalTangential_vs_layer",
      "Tangential component of track normal versus layer" + suffix +
      ";TPC layer;n_{r#phi}",
      48, 6.5, 54.5, 120, -1.2, 1.2);

    h.h_normalDistance_vs_layer = new TH2D(
      "h_normalDistance_vs_layer",
      "Signed normal displacement versus layer" + suffix +
      ";TPC layer;d_{normal} [cm]",
      48, 6.5, 54.5, 160, -maxAbsVoteCm, maxAbsVoteCm);

    h.h_closestPoint_deltaR_deltaRPhi = new TH2D(
      "h_closestPoint_deltaR_deltaRPhi",
      "Minimum-norm closest-point solution" + suffix +
      ";#Delta r [cm];r#Delta#phi [cm]",
      160, -maxAbsVoteCm, maxAbsVoteCm,
      160, -maxAbsVoteCm, maxAbsVoteCm);

    const int dimensions = 3;
    const int bins[dimensions] = {48, 100, 100};
    const double xmin[dimensions] = {6.5, -maxAbsVoteCm, -maxAbsVoteCm};
    const double xmax[dimensions] = {54.5, maxAbsVoteCm, maxAbsVoteCm};

    h.h_lineVotes_layer_deltaR_deltaRPhi = new THnSparseD(
      "h_lineVotes_layer_deltaR_deltaRPhi",
      "Allowed 2D correction-line votes;TPC layer;#Delta r [cm];r#Delta#phi [cm]",
      dimensions, bins, xmin, xmax);

    return h;
  }

  void fillAllowedCorrectionLine(
    THnSparseD* histogram,
    const double layer,
    const double normalR,
    const double normalPhi,
    const double normalDistance,
    const double maxAbsVoteCm,
    const int samples)
  {
    if (!histogram || samples < 2) return;

    // Fill the line n_r*DeltaR + n_phi*DeltaS = d.
    // Scan along the coordinate with the more stable denominator.
    if (std::abs(normalPhi) >= std::abs(normalR) && std::abs(normalPhi) > 1e-6)
    {
      for (int i = 0; i < samples; ++i)
      {
        const double deltaR = -maxAbsVoteCm +
          2. * maxAbsVoteCm * static_cast<double>(i) / (samples - 1);
        const double deltaS =
          (normalDistance - normalR * deltaR) / normalPhi;
        if (std::abs(deltaS) > maxAbsVoteCm) continue;
        const double values[3] = {layer, deltaR, deltaS};
        histogram->Fill(values, 1.0 / samples);
      }
    }
    else if (std::abs(normalR) > 1e-6)
    {
      for (int i = 0; i < samples; ++i)
      {
        const double deltaS = -maxAbsVoteCm +
          2. * maxAbsVoteCm * static_cast<double>(i) / (samples - 1);
        const double deltaR =
          (normalDistance - normalPhi * deltaS) / normalR;
        if (std::abs(deltaR) > maxAbsVoteCm) continue;
        const double values[3] = {layer, deltaR, deltaS};
        histogram->Fill(values, 1.0 / samples);
      }
    }
  }
}

void Make2DClusterCorrectionVotes(
  const char* inputDir = ".",
  const char* filePattern = "*.root",
  const char* outputName = "cluster_2d_correction_votes_v2.root",
  const char* treeName = "residuals",
  const double beamX = 0.158,
  const double beamY = 0.285,
  const double bFieldTesla = 1.4,
  const double r2AnchorMin = 46.0,
  const double r2AnchorMax = 54.0,
  const int nR2AnchorPoints = 9,
  const double maxAbsBeamDca = 2.0,
  const double maxAbsClusterZ = 15.0,
  const int nPhiBins = 12,
  const double maxAbsVoteCm = 0.8,
  const int lineSamples = 80,
  const Long64_t maxEntries = -1)
{
  TH1::AddDirectory(kTRUE);

  const TString chainPattern = TString::Format("%s/%s", inputDir, filePattern);
  TChain chain(treeName);
  const int nFiles = chain.Add(chainPattern);
  if (nFiles <= 0)
  {
    std::cerr << "ERROR: no files matched " << chainPattern << std::endl;
    return;
  }

  Int_t side = -1;
  UInt_t ntpc_clusters = 0;
  Double_t pt = 0., px = 0., py = 0., pz = 0.;
  Double_t charge = 0., quality = 0.;
  Double_t pca_x = 0., pca_y = 0., pca_z = 0.;

  std::vector<unsigned int>* layer = nullptr;
  std::vector<double>* cluster_r = nullptr;
  std::vector<double>* cluster_phi = nullptr;
  std::vector<double>* cluster_z = nullptr;

  chain.SetBranchAddress("side", &side);
  chain.SetBranchAddress("ntpc_clusters", &ntpc_clusters);
  chain.SetBranchAddress("pt", &pt);
  chain.SetBranchAddress("px", &px);
  chain.SetBranchAddress("py", &py);
  chain.SetBranchAddress("pz", &pz);
  chain.SetBranchAddress("charge", &charge);
  chain.SetBranchAddress("quality", &quality);
  chain.SetBranchAddress("pca_x", &pca_x);
  chain.SetBranchAddress("pca_y", &pca_y);
  chain.SetBranchAddress("pca_z", &pca_z);
  chain.SetBranchAddress("layer", &layer);
  chain.SetBranchAddress("cluster_r", &cluster_r);

  if (!chain.GetBranch("cluster_phi") || !chain.GetBranch("cluster_z"))
  {
    std::cerr << "ERROR: cluster_phi and cluster_z branches are required." << std::endl;
    return;
  }
  chain.SetBranchAddress("cluster_phi", &cluster_phi);
  chain.SetBranchAddress("cluster_z", &cluster_z);

  std::vector<double> anchorRadii;
  for (int i = 0; i < nR2AnchorPoints; ++i)
  {
    const double f = static_cast<double>(i) / (nR2AnchorPoints - 1);
    anchorRadii.push_back(r2AnchorMin + f * (r2AnchorMax - r2AnchorMin));
  }

  std::unique_ptr<TFile> output(TFile::Open(outputName, "RECREATE"));
  if (!output || output->IsZombie())
  {
    std::cerr << "ERROR: could not create " << outputName << std::endl;
    return;
  }

  std::map<std::string, GroupHistograms> groups;
  for (int iside = 0; iside < 2; ++iside)
  {
    for (int iq = 0; iq < 2; ++iq)
    {
      const double q = iq == 0 ? -1. : 1.;
      for (int ipt = 0; ipt < 4; ++ipt)
      {
        for (int iphi = 0; iphi < nPhiBins; ++iphi)
        {
          const std::string name = groupName(iside, q, ipt, iphi);
          TDirectory* dir = output.get();

          // Build nested directory path explicitly.
          TString path(name.c_str());
          TObjArray* tokens = path.Tokenize("/");
          for (int it = 0; it < tokens->GetEntries(); ++it)
          {
            const TString token = tokens->At(it)->GetName();
            TDirectory* next = dir->GetDirectory(token);
            if (!next) next = dir->mkdir(token);
            dir = next;
          }
          delete tokens;

          groups.emplace(name, bookGroup(dir, name, maxAbsVoteCm));
        }
      }
    }
  }

  output->cd();
  TH1D h_beamDcaBefore("h_beamDcaBefore",
    "Beam DCA before constraint;beam DCA [cm];tracks", 400, -10., 10.);
  TH1D h_r2AnchorRms("h_r2AnchorRms",
    "R2 anchor RMS;RMS [cm];tracks", 400, 0., 2.);
  TH1D h_absNormalRadial("h_absNormalRadial",
    "|radial component of track normal|;|n_{r}|;clusters", 100, 0., 1.);
  TH1D h_absNormalTangential("h_absNormalTangential",
    "|tangential component of track normal|;|n_{r#phi}|;clusters", 100, 0., 1.);

  const Long64_t totalEntries = chain.GetEntries();
  const Long64_t entriesToRun =
    maxEntries > 0 ? std::min(maxEntries, totalEntries) : totalEntries;

  Long64_t selectedTracks = 0;
  Long64_t selectedClusters = 0;

  for (Long64_t entry = 0; entry < entriesToRun; ++entry)
  {
    chain.GetEntry(entry);
    if (entry % 100000 == 0)
    {
      std::cout << "Processing " << entry << " / " << entriesToRun << std::endl;
    }

    const int ipt = ptBin(pt);
    if (ipt < 0 || (side != 0 && side != 1) ||
        !std::isfinite(pt) || !std::isfinite(px) || !std::isfinite(py) ||
        !std::isfinite(charge) || !std::isfinite(pca_x) ||
        !std::isfinite(pca_y) || !std::isfinite(pca_z) ||
        std::abs(charge) < 0.5 || ntpc_clusters <= 20 ||
        std::abs(pca_z) >= 10. || quality > 20.)
    {
      continue;
    }

    if (!layer || !cluster_r || !cluster_phi || !cluster_z) continue;

    const Circle original = buildTrackCircle(
      pca_x, pca_y, px, py, pt, charge, bFieldTesla);
    if (!original.valid) continue;

    const double beamDca =
      std::hypot(original.xc - beamX, original.yc - beamY) - original.radius;
    h_beamDcaBefore.Fill(beamDca);
    if (!std::isfinite(beamDca) || std::abs(beamDca) > maxAbsBeamDca) continue;

    const Circle constrained = constrainCircleToBeamAndR2(
      original, beamX, beamY, anchorRadii);
    if (!constrained.valid) continue;

    const double anchorCost = circleAnchorCost(original, constrained, anchorRadii);
    if (!std::isfinite(anchorCost)) continue;
    h_r2AnchorRms.Fill(std::sqrt(anchorCost));
    ++selectedTracks;

    const std::size_t nClusters = std::min(
      {layer->size(), cluster_r->size(), cluster_phi->size(), cluster_z->size()});

    for (std::size_t i = 0; i < nClusters; ++i)
    {
      const double r = cluster_r->at(i);
      const double phi = wrapPhi(cluster_phi->at(i));
      const double z = cluster_z->at(i);
      const double tpcLayer = static_cast<double>(layer->at(i));

      if (!std::isfinite(r) || !std::isfinite(phi) || !std::isfinite(z) ||
          r <= 0. || std::abs(z) >= maxAbsClusterZ)
      {
        continue;
      }

      const XYPoint measured = {r * std::cos(phi), r * std::sin(phi)};
      const XYPoint ideal = closestPointOnCircle(measured, constrained);
      if (!std::isfinite(ideal.x) || !std::isfinite(ideal.y)) continue;

      const double correctionX = ideal.x - measured.x;
      const double correctionY = ideal.y - measured.y;
      const double normalDistanceAbs = std::hypot(correctionX, correctionY);
      if (normalDistanceAbs < 1e-10 || normalDistanceAbs > maxAbsVoteCm) continue;

      // Local detector basis at the measured cluster.
      const double erX = std::cos(phi);
      const double erY = std::sin(phi);
      const double ephiX = -std::sin(phi);
      const double ephiY = std::cos(phi);

      const double deltaR = correctionX * erX + correctionY * erY;
      const double deltaRPhi = correctionX * ephiX + correctionY * ephiY;

      // Unit normal points from measured cluster toward the constrained track.
      const double normalX = correctionX / normalDistanceAbs;
      const double normalY = correctionY / normalDistanceAbs;
      const double normalR = normalX * erX + normalY * erY;
      const double normalPhi = normalX * ephiX + normalY * ephiY;
      const double normalDistance = normalDistanceAbs;

      h_absNormalRadial.Fill(std::abs(normalR));
      h_absNormalTangential.Fill(std::abs(normalPhi));

      const int iphi = phiBin(phi, nPhiBins);
      const std::string name = groupName(side, charge, ipt, iphi);
      GroupHistograms& h = groups.at(name);

      h.h_deltaR_vs_layer->Fill(tpcLayer, deltaR);
      h.h_deltaRPhi_vs_layer->Fill(tpcLayer, deltaRPhi);
      h.h_normalRadial_vs_layer->Fill(tpcLayer, normalR);
      h.h_normalTangential_vs_layer->Fill(tpcLayer, normalPhi);
      h.h_normalDistance_vs_layer->Fill(tpcLayer, normalDistance);
      h.h_closestPoint_deltaR_deltaRPhi->Fill(deltaR, deltaRPhi);

      fillAllowedCorrectionLine(
        h.h_lineVotes_layer_deltaR_deltaRPhi,
        tpcLayer,
        normalR,
        normalPhi,
        normalDistance,
        maxAbsVoteCm,
        lineSamples);

      ++selectedClusters;
    }
  }

  output->cd();
  output->Write();
  output->Close();

  std::cout << "Selected tracks: " << selectedTracks << '\n'
            << "Selected |z_cluster|<" << maxAbsClusterZ
            << " cm clusters: " << selectedClusters << '\n'
            << "Wrote: " << outputName << std::endl;
}
