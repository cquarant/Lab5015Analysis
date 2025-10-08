import ROOT
import sys


# List of ROOT files to process
# You can also pass them as arguments: python average_graph_y.py file1.root file2.root ...
file_dict = [
    { 
        "DM" : "FE4587",
        "file_list": ["/home/cmsdaq/TBAnalysis/TB_CERN_Sept2025/Lab5015Analysis/plots/summaryPlots_ith1_scan_FE4587_run4269.root",
                      "/home/cmsdaq/TBAnalysis/TB_CERN_Sept2025/Lab5015Analysis/plots/summaryPlots_ith1_scan_FE4587_run4272.root",
                      "/home/cmsdaq/TBAnalysis/TB_CERN_Sept2025/Lab5015Analysis/plots/summaryPlots_analysis_std_run4278.root",
                      "/home/cmsdaq/TBAnalysis/TB_CERN_Sept2025/Lab5015Analysis/plots/summaryPlots_analysis_std_run4281.root"], 
        "ovs": ["0.90", "1.20", "2.00", "3.00"],
        "Vbd":38.04, 
        "color":ROOT.kBlack,
        "marker":ROOT.kFullCircle,
    },
    { 
        "DM": "DM9000",
        "file_list": ["/home/cmsdaq/TBAnalysis/TB_CERN_Sept2025/Lab5015Analysis/plots/summaryPlots_ith1_scan_FE3966_run4140.root",
                      "/home/cmsdaq/TBAnalysis/TB_CERN_Sept2025/Lab5015Analysis/plots/summaryPlots_ith1_scan_FE3966_run4141-4143.root",
                      "/home/cmsdaq/TBAnalysis/TB_CERN_Sept2025/Lab5015Analysis/plots/summaryPlots_ith1_scan_FE3966_run4150.root",
                      "/home/cmsdaq/TBAnalysis/TB_CERN_Sept2025/Lab5015Analysis/plots/summaryPlots_ith1_scan_FE3966_run4166.root"],
        "ovs": ["3.00", "2.00", "1.20", "0.90"],
        "Vbd":38.56,
        "color":ROOT.kRed,
        "marker":ROOT.kFullSquare,
    }
]

setVbd = 38.11
tRes_vs_OV_dict = {}
tGraphs = {}
for d in file_dict:

    tRes_vs_OV_dict[d['DM']] = {'tRes':[], 'ovs':[]}
    tgraph = ROOT.TGraph()
    tgraph.SetName("tRes_vs_OV_"+d['DM'])
    
    print d['DM']
    for j,filename in enumerate(d['file_list']):
        
        graph_name = "g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_bar_Vov"+str(d['ovs'][j])+"_enBin01"

        f = ROOT.TFile.Open(filename)
        if not f or f.IsZombie():
            print "Could not open "+filename
            continue

        g = f.Get(graph_name)
        if not g:
            print "Graph '{}' not found in {}".format(graph_name, filename)
            f.Close()
            continue

        n = g.GetN()
        if n == 0:
            print "Graph '{}' has no points in {}".format(graph_name, filename)
            f.Close()
            continue

        # Compute mean of Y values
        y_values = [g.GetPointY(i) for i in range(n)]
        mean_y = sum(y_values) / n

        tRes_vs_OV_dict[d['DM']]['tRes'].append(mean_y)
        tRes_vs_OV_dict[d['DM']]['ovs'].append(float(d['ovs'][j]) + (d['Vbd'] - setVbd))
        tgraph.SetPoint(j, float(d['ovs'][j]) - (d['Vbd'] - setVbd), mean_y)
        
        print tgraph.GetPointY(j), " at OV ", tgraph.GetPointX(j)


    tGraphs[d['DM']] = tgraph

c = ROOT.TCanvas("c", "c", 800, 600)
c.cd()

#Dummy histogram for axes
h_dummy = ROOT.TH2F("h_dummy", "", 10, 0, 3.5, 10, 0, 100)
h_dummy.Draw()

for i, gr in enumerate(tGraphs.keys()):
    tGraphs[gr].SetMarkerStyle(file_dict[i]['marker'])
    tGraphs[gr].SetMarkerColor(file_dict[i]['color'])

    tGraphs[gr].Draw("P SAME")

c.Update()
c.SaveAs("timeRes_vs_OVcorrected.png")
c.SaveAs("timeRes_vs_OVcorrected.pdf")