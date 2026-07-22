// Make2DClusterCorrectionVotes_v3.C
//
// Test of a common two-dimensional transverse cluster correction.
//
// Selection and directory structure:
//   side / charge / pT / module / phiLayerBin
//
// pT bins:
//   0.2-0.4, 0.4-0.7, 0.7-1.2, 1.2-1.8, 1.8-5.0 GeV/c
//
// A hardware module is identified by:
//   12 azimuthal sectors x 3 radial modules = 36 modules per side.
// Each module is divided into:
//   3 local-phi bins x 5 layer bins.
// The first and last layer of each radial module are excluded.
//
// Only one TH2D is booked in every final bin:
//   x = Delta r [cm]
//   y = r Delta phi [cm]
//
// IMPORTANT:
// A cluster does not fill its closest point on the constrained track. Instead,
// it fills the complete allowed correction line
//
//   n_r Delta r + n_phi (r Delta phi) = d_normal.
//
// Each cluster contributes the same total weight along the visible part of its
// allowed line. Therefore high-pT tracks can remain nearly uniform in Delta r,
// while lower-pT tracks provide tilted lines. A common 2D maximum can appear
// only from the intersection of many such lines.
//
// Example:
// root -l -b -q 'Make2DClusterCorrectionVotes_v3.C("/path","*.root","votes_v3.root")'

#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TString.h>

#include <algorithm>
#include <array>
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

  double wrapPhi(double phi)
  {
    while (phi >= TMath::Pi()) phi -= kTwoPi;
    while (phi < -TMath::Pi()) phi += kTwoPi;
    return phi;
  }

  double wrapToTwoPi(double phi)
  {
    phi = std::fmod(phi, kTwoPi);
    if (phi < 0.) phi += kTwoPi;
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

    for (const double radius : anchorRadii)
    {
      const auto oldPoints = intersectCircleWithOriginRadius(original, radius);
      const auto newPoints = intersectCircleWithOriginRadius(candidate, radius);
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
      constexpr int nLocal = 80;
      double localBestTheta = bestTheta;
      double localBestCost = bestCost;

      for (int i = 0; i <= nLocal; ++i)
      {
        const double fraction = static_cast<double>(i) / nLocal;
        const double theta = bestTheta - halfWidth + 2. * halfWidth * fraction;

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
    if (pt >= 0.2 && pt < 0.23) return 0;
    if (pt >= 0.30 && pt < 0.4) return 1;
    if (pt >= 0.4 && pt < 0.45) return 2;
    if (pt >= 0.45 && pt < 0.7) return 3;
    if (pt >= 0.7 && pt < 5.0) return 4;
    return -1;
  }

  std::string ptBinName(const int bin)
  {
    static const std::array<const char*, 5> names = {
      "pt_0p2_0p4",
      "pt_0p4_0p7",
      "pt_0p7_1p2",
      "pt_1p2_1p8",
      "pt_1p8_5p0"
    };
    return (bin >= 0 && bin < static_cast<int>(names.size()))
      ? names[bin] : "pt_invalid";
  }

  // Radial module layer ranges: R1=7-22, R2=23-38, R3=39-54.
  int radialModuleFromLayer(const int layer)
  {
    if (layer >= 7 && layer <= 22) return 0;
    if (layer >= 23 && layer <= 38) return 1;
    if (layer >= 39 && layer <= 54) return 2;
    return -1;
  }

  int firstLayerForRadialModule(const int radialModule)
  {
    return 7 + 16 * radialModule;
  }

  // Exclude the first and last layer. The remaining 14 layers are divided
  // as evenly as possible into five bins: 3,3,3,3,2 layers.
  int localLayerBin(const int layer, const int radialModule)
  {
    if (radialModule < 0 || radialModule > 2) return -1;
    const int first = firstLayerForRadialModule(radialModule);
    const int local = layer - first;
    if (local <= 0 || local >= 15) return -1;

    const int index = local - 1; // 0..13 after removing edge layers
    if (index < 3) return 0;
    if (index < 6) return 1;
    if (index < 9) return 2;
    if (index < 12) return 3;
    return 4;
  }

  struct PhiModuleBin
  {
    int sector = -1;
    int localPhiBin = -1;
    double localFraction = 0.;
  };

  // sectorPhiOffset is the low edge of sector 0. The default -pi gives 12
  // uniform sectors over [-pi,pi). Adjust it if the hardware sector boundaries
  // use another phase convention.
  PhiModuleBin determinePhiModuleBin(
    const double phi,
    const double sectorPhiOffset)
  {
    PhiModuleBin result;
    const double sectorWidth = kTwoPi / 12.;
    const double shifted = wrapToTwoPi(phi - sectorPhiOffset);
    int sector = static_cast<int>(shifted / sectorWidth);
    sector = std::clamp(sector, 0, 11);

    const double withinSector = shifted - sector * sectorWidth;
    const double fraction = withinSector / sectorWidth;
    int localPhi = static_cast<int>(3. * fraction);
    localPhi = std::clamp(localPhi, 0, 2);

    result.sector = sector;
    result.localPhiBin = localPhi;
    result.localFraction = fraction;
    return result;
  }

  std::string groupName(
    const int side,
    const double charge,
    const int iPt,
    const int hardwareModule,
    const int localPhi,
    const int localLayer)
  {
    return TString::Format(
      "side%d/%s/%s/module_%02d/phi_%d_layer_%d",
      side,
      charge > 0. ? "qplus" : "qminus",
      ptBinName(iPt).c_str(),
      hardwareModule,
      localPhi,
      localLayer).Data();
  }

  TDirectory* makeDirectoryPath(TDirectory* root, const std::string& path)
  {
    TDirectory* current = root;
    TString temporary(path.c_str());
    std::unique_ptr<TObjArray> tokens(temporary.Tokenize("/"));

    for (int i = 0; i < tokens->GetEntries(); ++i)
    {
      const TString token = tokens->At(i)->GetName();
      TDirectory* next = current->GetDirectory(token);
      if (!next) next = current->mkdir(token);
      current = next;
    }
    return current;
  }

  TH2D* bookVoteHistogram(
    TDirectory* directory,
    const std::string& label,
    const int voteBins,
    const double maxAbsVoteCm)
  {
    directory->cd();
    TH2D* histogram = new TH2D(
      "h_deltaRPhi_vs_deltaR_votes",
      TString::Format(
        "Allowed 2D correction votes [%s];#Delta r [cm];r#Delta#phi [cm]",
        label.c_str()),
      voteBins, -maxAbsVoteCm, maxAbsVoteCm,
      voteBins, -maxAbsVoteCm, maxAbsVoteCm);
    histogram->Sumw2();
    return histogram;
  }

  // Fill the visible part of one allowed line with unit total weight.
  // This avoids a fake preference for Delta r = 0 and avoids giving more
  // weight to lines that happen to cross more histogram bins.
  void fillAllowedCorrectionLine(
    TH2D* histogram,
    const double normalR,
    const double normalPhi,
    const double normalDistance,
    const double maxAbsVoteCm,
    const int lineSamples)
  {
    if (!histogram || lineSamples < 2) return;

    std::vector<std::pair<double, double>> points;
    points.reserve(lineSamples);

    if (std::abs(normalPhi) >= std::abs(normalR) && std::abs(normalPhi) > 1e-8)
    {
      for (int i = 0; i < lineSamples; ++i)
      {
        const double deltaR = -maxAbsVoteCm +
          2. * maxAbsVoteCm * static_cast<double>(i) / (lineSamples - 1);
        const double deltaRPhi =
          (normalDistance - normalR * deltaR) / normalPhi;
        if (std::abs(deltaRPhi) <= maxAbsVoteCm)
        {
          points.emplace_back(deltaR, deltaRPhi);
        }
      }
    }
    else if (std::abs(normalR) > 1e-8)
    {
      for (int i = 0; i < lineSamples; ++i)
      {
        const double deltaRPhi = -maxAbsVoteCm +
          2. * maxAbsVoteCm * static_cast<double>(i) / (lineSamples - 1);
        const double deltaR =
          (normalDistance - normalPhi * deltaRPhi) / normalR;
        if (std::abs(deltaR) <= maxAbsVoteCm)
        {
          points.emplace_back(deltaR, deltaRPhi);
        }
      }
    }

    if (points.empty()) return;
    const double weight = 1.0 / static_cast<double>(points.size());
    for (const auto& point : points)
    {
      histogram->Fill(point.first, point.second, weight);
    }
  }
}

void Make2DClusterCorrectionVotes(
  const char* inputDir = ".",
  const char* filePattern = "*.root",
  const char* outputName = "cluster_2d_correction_votes_v3.root",
  const char* treeName = "residuals",
  const double beamX = 0.158,
  const double beamY = 0.285,
  const double bFieldTesla = 1.4,
  const double r2AnchorMin = 46.0,
  const double r2AnchorMax = 54.0,
  const int nR2AnchorPoints = 9,
  const double maxAbsBeamDca = 2.0,
  const double maxAbsClusterZ = 15.0,
  const double sectorPhiOffset = -TMath::Pi(),
  const double maxAbsVoteCm = 0.8,
  const int voteBins = 121,
  const int lineSamples = 241,
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
  if (nR2AnchorPoints < 2)
  {
    std::cerr << "ERROR: nR2AnchorPoints must be at least 2." << std::endl;
    return;
  }
  for (int i = 0; i < nR2AnchorPoints; ++i)
  {
    const double fraction = static_cast<double>(i) / (nR2AnchorPoints - 1);
    anchorRadii.push_back(r2AnchorMin + fraction * (r2AnchorMax - r2AnchorMin));
  }

  std::unique_ptr<TFile> output(TFile::Open(outputName, "RECREATE"));
  if (!output || output->IsZombie())
  {
    std::cerr << "ERROR: could not create " << outputName << std::endl;
    return;
  }

  std::map<std::string, TH2D*> voteHistograms;

  output->cd();
  TH1D h_beamDcaBefore(
    "h_beamDcaBefore",
    "Beam DCA before constraint;beam DCA [cm];tracks",
    400, -10., 10.);
  TH1D h_r2AnchorRms(
    "h_r2AnchorRms",
    "R2 anchor RMS;RMS [cm];tracks",
    400, 0., 2.);
  TH1D h_absNormalRadial(
    "h_absNormalRadial",
    "|radial component of track normal|;|n_{r}|;clusters",
    100, 0., 1.);
  TH1D h_absNormalTangential(
    "h_absNormalTangential",
    "|tangential component of track normal|;|n_{r#phi}|;clusters",
    100, 0., 1.);
  TH1D h_clustersPerVoteHistogram(
    "h_clustersPerVoteHistogram",
    "Number of cluster votes per final histogram;cluster votes;histograms",
    200, 0., 20000.);

  std::map<std::string, Long64_t> voteCounts;

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

    const int iPt = ptBin(pt);
    if (iPt < 0 || (side != 0 && side != 1) ||
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
      const int tpcLayer = static_cast<int>(layer->at(i));
      const int radialModule = radialModuleFromLayer(tpcLayer);
      const int iLayer = localLayerBin(tpcLayer, radialModule);
      if (radialModule < 0 || iLayer < 0) continue;

      const double radius = cluster_r->at(i);
      const double phi = wrapPhi(cluster_phi->at(i));
      const double z = cluster_z->at(i);

      if (!std::isfinite(radius) || !std::isfinite(phi) || !std::isfinite(z) ||
          radius <= 0. || std::abs(z) >= maxAbsClusterZ)
      {
        continue;
      }

      const PhiModuleBin phiInfo = determinePhiModuleBin(phi, sectorPhiOffset);
      if (phiInfo.sector < 0 || phiInfo.localPhiBin < 0) continue;

      // 36 hardware modules per side: three radial modules in each sector.
      const int hardwareModule = 3 * phiInfo.sector + radialModule;

      const XYPoint measured = {
        radius * std::cos(phi),
        radius * std::sin(phi)
      };
      const XYPoint closest = closestPointOnCircle(measured, constrained);
      if (!std::isfinite(closest.x) || !std::isfinite(closest.y)) continue;

      const double correctionX = closest.x - measured.x;
      const double correctionY = closest.y - measured.y;
      const double normalDistance = std::hypot(correctionX, correctionY);
      if (!std::isfinite(normalDistance) || normalDistance < 1e-10 ||
          normalDistance > maxAbsVoteCm)
      {
        continue;
      }

      const double radialX = std::cos(phi);
      const double radialY = std::sin(phi);
      const double tangentialX = -std::sin(phi);
      const double tangentialY = std::cos(phi);

      // Unit normal points from the measured cluster to the constrained track.
      const double normalX = correctionX / normalDistance;
      const double normalY = correctionY / normalDistance;
      const double normalR = normalX * radialX + normalY * radialY;
      const double normalPhi = normalX * tangentialX + normalY * tangentialY;

      h_absNormalRadial.Fill(std::abs(normalR));
      h_absNormalTangential.Fill(std::abs(normalPhi));

      const std::string name = groupName(
        side, charge, iPt, hardwareModule,
        phiInfo.localPhiBin, iLayer);

      TH2D*& histogram = voteHistograms[name];
      if (!histogram)
      {
        TDirectory* directory = makeDirectoryPath(output.get(), name);
        histogram = bookVoteHistogram(
          directory, name, voteBins, maxAbsVoteCm);
      }

      fillAllowedCorrectionLine(
        histogram,
        normalR,
        normalPhi,
        normalDistance,
        maxAbsVoteCm,
        lineSamples);

      ++voteCounts[name];
      ++selectedClusters;
    }
  }

  for (const auto& [name, count] : voteCounts)
  {
    h_clustersPerVoteHistogram.Fill(static_cast<double>(count));
  }

  output->cd();
  output->Write();
  output->Close();

  std::cout << "Selected tracks: " << selectedTracks << '\n'
            << "Selected clusters: " << selectedClusters << '\n'
            << "Booked final vote histograms: " << voteHistograms.size() << '\n'
            << "Wrote: " << outputName << std::endl;
}
