from __future__ import print_function

__author__ = "Rita Giordano, Gianluca Santoni, Tiankun Zhou"
__copyright__ = "Copyright 2015-2026"
__credits__ = ["Rita Giordano, Gianluca Santoni, Tiankun Zhou, Alexander Popov"]
__license__ = ""
__version__ = "2.0"
__maintainer__ = "Tiankun Zhou"
__email__ = "tiankun.zhou@esrf.fr"
__status__ = "Beta"


from PySide6 import QtCore, QtWidgets
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.lines import Line2D
from pathlib import Path
from scipy.cluster import hierarchy
from .clustering import extractXSCALEStat, checkIndices, find_cluster_for_index, colors
import collections
import operator

# implement the default mpl key bindings
import matplotlib.patches as mpatches
import numpy as np

#create results tab to plot the SCALE.LP statistics
class SinglePlotTab(QtWidgets.QWidget):
    def __init__(self, ProcessedDir:str):
        super().__init__()

        self.ProcessedDir= ProcessedDir
        #setup widget
        self.tabLayout= QtWidgets.QVBoxLayout(self)
        self.Title=QtWidgets.QLabel(self)
        ResultName = Path(self.ProcessedDir).name.split('_')
        self.Title.setText(f"Threshold: {ResultName[2]} and Group number: {ResultName[3]}")
        self.Title.setStyleSheet("color: black; font-weight: bold; font-size: 14px;")

        #set up the plot
        self.statsPlot = Figure()
        self.Ax = self.statsPlot.add_subplot(111)
        self.statsCanvas = FigureCanvas(self.statsPlot)
        self.statsBar= NavigationToolbar(self.statsCanvas, self)

        #Buttons to plot stats
        self.buttonBar = QtWidgets.QWidget()
        self.barLayout= QtWidgets.QHBoxLayout(self.buttonBar)

        #plot from XSCLE.LP, position 0 is resolution, 4 is COMPLETENESS
        #5 is R factor observed, 8 is <I/sigma>, 10 is CC1/2, 12 is SigAno
        self.ccVsR= QtWidgets.QPushButton()
        self.ccVsR.setText("CC. vs Res")
        self.ccVsR.clicked.connect(lambda:self.plotStats(0, 10, "CC. vs Res"))

        self.compVsR= QtWidgets.QPushButton()
        self.compVsR.setText("comp vs Res")
        self.compVsR.clicked.connect(lambda:self.plotStats(0, 4, "comp vs Res"))

        self.RobsVsR= QtWidgets.QPushButton()
        self.RobsVsR.setText("Robs vs Res")
        self.RobsVsR.clicked.connect(lambda:self.plotStats(0, 5, "Robs vs Res"))

        self.IsigmaVsR= QtWidgets.QPushButton()
        self.IsigmaVsR.setText("<I/\u03C3I> vs Res")
        self.IsigmaVsR.clicked.connect(lambda:self.plotStats(0, 8, "<I/\u03C3I> vs Res"))

        self.SanoVsR= QtWidgets.QPushButton()
        self.SanoVsR.setText("Sig. Ano. vs Res")
        self.SanoVsR.clicked.connect(lambda:self.plotStats(0, 12, "Sig. Ano. vs Res"))

        self.barLayout.addWidget(self.ccVsR)
        self.barLayout.addWidget(self.compVsR)
        self.barLayout.addWidget(self.RobsVsR)
        self.barLayout.addWidget(self.SanoVsR)
        self.barLayout.addWidget(self.IsigmaVsR)

        self.tabLayout.addWidget(self.Title, 1, alignment=QtCore.Qt.AlignCenter)
        self.tabLayout.addWidget(self.buttonBar, 1)
        self.tabLayout.addWidget(self.statsBar, 1)
        self.tabLayout.addWidget(self.statsCanvas, 17)
        
    def plotStats(self, res, value, title):
        self.Ax.clear()
        plotDataX= []
        plotDataY= []
        plotList, _ = extractXSCALEStat(f"{self.ProcessedDir}/XSCALE.LP")

        #Check if plotList is empty
        if not plotList:
            self.Ax.text(0.5, 0.5, "No data available to plot", 
                         horizontalalignment='center', verticalalignment='center')
            self.statsCanvas.draw()
            return

        LowestRes = float(plotList[0][res]) + 0.25
        HighestRes = float(plotList[-1][res]) - 0.25 if float(plotList[-1][res]) > 0.25 else 0
        for line in plotList:
            plotDataX.append(float(line[res])) 
            plotDataY.append(float(line[value].strip('*').strip('%')))
        self.Ax.plot(plotDataX, plotDataY, 'r-^')
        self.Ax.set_xlim(LowestRes, HighestRes)
        #X and y labels
        self.Ax.set_xlabel("Resolution (Å)", fontsize=12, fontweight='bold', color='black')
        y_label = title.split(" vs ")[0]
        self.Ax.set_ylabel(y_label, fontsize=12, fontweight='bold', rotation=90, labelpad=20, color='black')
        self.Ax.set_title(title)

        self.statsPlot.tight_layout()
        self.statsCanvas.draw()


#plot multple statistics for multiple processed directories
class MultiPlotTab(QtWidgets.QWidget):
    def __init__(self, ProcessedDirs:list):
        super().__init__()
        self.ProcessedDirs= ProcessedDirs
        #set up widget
        self.tabLayout= QtWidgets.QVBoxLayout(self)
        self.Title=QtWidgets.QLabel(self)
        self.Title.setText(f"Multi statistics for {len(self.ProcessedDirs)} processed directories")
        self.Title.setStyleSheet("color: black; font-weight: bold; font-size: 14px;")

        #set up the plot
        self.MultiStatusPlots = Figure()
        self.Ax = self.MultiStatusPlots.add_subplot(111)
        self.MultiStatsCanvas = FigureCanvas(self.MultiStatusPlots)
        self.MultiplestatsBar= NavigationToolbar(self.MultiStatsCanvas, self)

        #Buttons to plot stats
        self.buttonBar = QtWidgets.QWidget()
        self.barLayout= QtWidgets.QHBoxLayout(self.buttonBar)

        self.ccVsR= QtWidgets.QPushButton()
        self.ccVsR.setText("CC. vs Res")
        self.ccVsR.clicked.connect(lambda:self.MultiPlotStats(0, 10, "CC. vs Res Multi"))

        self.compVsR= QtWidgets.QPushButton()
        self.compVsR.setText("comp vs Res")
        self.compVsR.clicked.connect(lambda:self.MultiPlotStats(0, 4, "comp vs Res Multi"))

        self.RobsVsR= QtWidgets.QPushButton()
        self.RobsVsR.setText("Robs vs Res")
        self.RobsVsR.clicked.connect(lambda:self.MultiPlotStats(0, 5, "Robs vs Res Multi"))

        self.IsigmaVsR= QtWidgets.QPushButton()
        self.IsigmaVsR.setText("<I/\u03C3I> vs Res")
        self.IsigmaVsR.clicked.connect(lambda:self.MultiPlotStats(0, 8, "<I/\u03C3I> vs Res Multi"))

        self.SanoVsR= QtWidgets.QPushButton()
        self.SanoVsR.setText("Sig. Ano. vs Res")
        self.SanoVsR.clicked.connect(lambda:self.MultiPlotStats(0, 12,"Sig. Ano. vs Res Multi"))

        self.barLayout.addWidget(self.ccVsR)
        self.barLayout.addWidget(self.compVsR)
        self.barLayout.addWidget(self.RobsVsR)
        self.barLayout.addWidget(self.SanoVsR)
        self.barLayout.addWidget(self.IsigmaVsR)

        self.tabLayout.addWidget(self.Title, 1, alignment=QtCore.Qt.AlignCenter)
        self.tabLayout.addWidget(self.buttonBar, 1)
        self.tabLayout.addWidget(self.MultiplestatsBar, 1)
        self.tabLayout.addWidget(self.MultiStatsCanvas, 17)


    #plot the statistics for the selected button
    def MultiPlotStats(self, res, value, title):
        self.Ax.clear()

        #store legends
        legend_handles = []

        #setup plot for each processed directory
        for ProcessedDir in self.ProcessedDirs:
            plotList, _ = extractXSCALEStat(f"{ProcessedDir}/XSCALE.LP")

            #Check if plotList is empty
            if not plotList:
                self.Ax.text(0.5, 0.5, f"No data available to plot for {Path(ProcessedDir).name}", 
                             horizontalalignment='center', verticalalignment='center')
                continue

            plotDataX= []
            plotDataY= []
            LowestResList = []
            HighestResList = []
            
            for line in plotList:
                #check whether x and v values are valid numbers
                try:
                    LowestResList.append(float(line[res]) + 0.25)
                    HighestResList.append(float(line[res]) - 0.25 if float(line[res]) > 0.5 else 0)
                    x_value = float(line[res]) if line[res] else 0
                    y_value = float(line[value].strip('*').strip('%')) if line[value] else 0
                    plotDataX.append(x_value)
                    plotDataY.append(y_value)
                except (ValueError, TypeError) as e:
                    print(f"Error converting data: {e} - line[{value}] = {line[value]}")
                    continue

            #add figure legend with the name of the processed directory
            line, = self.Ax.plot(plotDataX, plotDataY, '-^', label=Path(ProcessedDir).name)
            legend_handles.append(line)

        #X and y labels
        self.Ax.set_xlabel("Resolution (Å)", fontsize=12, fontweight='bold', color='black')
        y_label = title.split(" vs ")[0]
        self.Ax.set_ylabel(y_label, fontsize=12, fontweight='bold', rotation=90, labelpad=20, color='black')

        self.Ax.set_xlim(max(LowestResList), min(HighestResList))
        self.Ax.set_title(title)

        #setup the legend
        if legend_handles:
            self.Ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.02, 1),\
                           borderaxespad=0, fontsize="small", title_fontsize="medium")
            self.MultiStatusPlots.tight_layout()
        else:
            self.MultiStatusPlots.tight_layout()

        #plot the figure
        self.MultiStatsCanvas.draw()



#create pre-plot dendrogram tab to plot the dendrogram with given Threshold
class PrePlotDendrogram(QtWidgets.QWidget):
    def __init__(self, ccClusterfile:QtWidgets.QLineEdit, Threshold:QtWidgets.QLineEdit, SelectedCluster:QtWidgets.QLineEdit, setupCC_method):
        super().__init__()
        self.ccClusterfile_widget = ccClusterfile
        self.threshold_widget = Threshold
        self.setupCC_method = setupCC_method
        self.SelectedCluster = SelectedCluster

        #setup layout
        self.tabLayout = QtWidgets.QVBoxLayout(self)
        self.tabLayout.setContentsMargins(0, 0, 0, 0)

        #setup tabs widget
        self.dendroClusterTabs = QtWidgets.QTabWidget()

        #setup dendro tabs
        self.dendrotab = QtWidgets.QWidget()
        self.dendrotabLayout = QtWidgets.QVBoxLayout(self.dendrotab)

        self.dendroTitle=QtWidgets.QLabel(self)
        self.dendroTitle.setText(f"<html><span style='color: black; font-weight: bold; font-size: 14px;'>Dendrogram with threshold: \
                                 {self.threshold_widget.text().strip() if self.threshold_widget else ''}</span></html>")
    
        #set up the plot
        self.dendroPlot = Figure()
        self.Ax = self.dendroPlot.add_subplot(111)
        self.dendroCanvas = FigureCanvas(self.dendroPlot)
        self.dendrostatsBar= NavigationToolbar(self.dendroCanvas, self)

        #plot button
        self.dendroplotButton = QtWidgets.QPushButton("Plot Dendrogram")
        self.dendroplotButton.setFixedSize(200, 30)
        self.dendroplotButton.clicked.connect(self.on_plot_clicked)

        #add Dendrogram widgets to layout
        self.dendrotabLayout.addWidget(self.dendroTitle, 1, alignment=QtCore.Qt.AlignCenter)
        self.dendrotabLayout.addWidget(self.dendroplotButton, 1, alignment=QtCore.Qt.AlignCenter)
        self.dendrotabLayout.addWidget(self.dendrostatsBar, 1)
        self.dendrotabLayout.addWidget(self.dendroCanvas, 27)

        #Set up the clutster content tab
        self.clusterContentTab = QtWidgets.QWidget()
        self.clusterContentLayout = QtWidgets.QVBoxLayout(self.clusterContentTab)

        #setup the cluster content text area
        self.clustercontentTitle=QtWidgets.QLabel(self)
        self.clustercontentTitle.setText(f"<html><span style='color: black; font-weight: bold; font-size: 14px;'>Cluster content with threshold: {self.threshold_widget.text().strip() if self.threshold_widget else ''}\
                                        The merged cluster and the corresponding datasets will be shown in <span style='color: blue;'>BLUE</span> color</span></html>")
        self.ShowNUmAndPathTitle=QtWidgets.QLabel(self)
        self.ShowNUmAndPathTitle.setText(f"<html><span style='color: black; font-weight: bold; font-size: 12px;'>The number and path of the datasets in \
                                        the merged cluster will be shown below as HTML for checking</span></html>")
        self.clusterContentText = QtWidgets.QTextEdit()
        self.clusterContentText.setReadOnly(True)
        self.clusterContentText.setStyleSheet("background-color: #f0f0f0; font-size: 14px;")
        self.ShowOnlyPathTitle=QtWidgets.QLabel(self)
        self.ShowOnlyPathTitle.setText(f"<html><span style='color: black; font-weight: bold; font-size: 12px;'>The path of the datasets in the \
                                      merged cluster will be shown below as plain text for copying</span></html>")
        self.clusterPathText = QtWidgets.QPlainTextEdit()
        self.clusterPathText.setReadOnly(True)
        self.clusterPathText.setStyleSheet("background-color: #f0f0f0; font-size: 14px;")

        #add the cluster content text area to the layout
        self.clusterContentLayout.addWidget(self.clustercontentTitle, 1, alignment=QtCore.Qt.AlignCenter)
        self.clusterContentLayout.addWidget(self.ShowNUmAndPathTitle, 1, alignment=QtCore.Qt.AlignCenter)
        self.clusterContentLayout.addWidget(self.clusterContentText, 28)
        self.clusterContentLayout.addWidget(self.ShowOnlyPathTitle, 1, alignment=QtCore.Qt.AlignCenter)
        self.clusterContentLayout.addWidget(self.clusterPathText, 14)

        #add tabs 
        self.dendroClusterTabs.addTab(self.dendrotab, "Dendrogram Plot")
        self.dendroClusterTabs.addTab(self.clusterContentTab, "Cluster Content")

        #add widget
        self.tabLayout.addWidget(self.dendroClusterTabs)


    #check whether the threshold is valid float and plot the dendrogram
    def on_plot_clicked(self):
        #get current values from the widgets
        ccClusterfile = self.ccClusterfile_widget.text() if self.ccClusterfile_widget else None
        threshold_text = self.threshold_widget.text().strip() if self.threshold_widget else None
        SelectedClusterText = self.SelectedCluster.text().replace(',', ' ').strip() if self.SelectedCluster else None
        
        if threshold_text:
            try:
                threshold = round(float(threshold_text), 2)
                self.PlotDendrogram(ccClusterfile, threshold)
                self.showMergedCluster(ccClusterfile, threshold, SelectedClusterText)
            except ValueError:
                #Show error for invalid threshold
                self.Ax.clear()
                self.Ax.text(0.5, 0.5, "Invalid threshold value", 
                            horizontalalignment='center', verticalalignment='center')
                self.dendroCanvas.draw()
        else:
            self.PlotDendrogram(ccClusterfile, None)
            self.showMergedCluster(ccClusterfile, None, SelectedClusterText)


    #prepare Dendrogram with current threshold
    def PlotDendrogram(self, ccClusterfile:str=None, threshold:float=None):
        #Clean plot
        self.Ax.clear()
        
        # Check parameters
        if ccClusterfile is None:
            self.Ax.text(0.5, 0.5, "No ccCluster file in Work dir", \
                         horizontalalignment='center', verticalalignment='center')
            self.dendroCanvas.draw()
            return None
        
        if threshold is None:
            self.Ax.text(0.5, 0.5, "No threshold provided", \
                         horizontalalignment='center', verticalalignment='center')
            self.dendroCanvas.draw()
            return None

        #set up CC
        _, Tree, etiquets, _ = self.setupCC_method(ccClusterfile)

        #check whether the dendrogram tree is valid
        if Tree is None:
            self.Ax.text(0.5, 0.5, "Error: Could not generate dendrogram",
                        horizontalalignment='center', verticalalignment='center')
            self.dendroCanvas.draw()
            return None

        X = hierarchy.dendrogram(Tree, color_threshold=threshold, ax=self.Ax)

        #try to match the color to cluster number using count
        #As they contains same amount od datasets ONLY IF the number of cluster is smaller than 10
        #get the FlatC with the same threshold
        FlatC = hierarchy.fcluster(Tree, threshold, criterion='distance')

        #Count the occurrences of each cluster and color
        cluster_counts = collections.Counter(FlatC)
        color_counts = collections.Counter(X['leaves_color_list'])

        #Remove unwanted color
        for remove_color in ['C0', 'k', 'grey']:
            if remove_color in color_counts:
                del color_counts[remove_color]

        #Remove cluster with only 1 dataset
        filtered_clusters = {cluster: count for cluster, count in cluster_counts.items() if count > 1}

        #Sort by count from largest to smallest, prepare to show maximum FIVE largest clusters in the legend
        sorted_clusters = sorted(filtered_clusters.items(), key=lambda x: x[1], reverse=True)
        sorted_colors = sorted(color_counts.items(), key=lambda x: x[1], reverse=True)
        print(f"sorted cluster: {sorted_clusters}")
        print(f"sorted color: {sorted_colors}")

        #Match clusters to colors by count
        cluster_to_color = {}
        for (cluster, cluster_count), (color, color_count) in zip(sorted_clusters, sorted_colors):
            if cluster_count == color_count:
                cluster_to_color[cluster] = color
            else:
                print(f"Warning: Count mismatch! Cluster {cluster} has {cluster_count}, Color {color} has {color_count}, possible due to large number of clusters")
                cluster_to_color[cluster] = color

        #Create legend using matched clusters
        legend_handles = []

        #add x and y labels, as well as the threshold line
        self.Ax.set_xlabel("Dateset numbers", fontsize=12, fontweight='bold', color='black')
        self.Ax.set_ylabel("Correlation coefficients", fontsize=12, fontweight='bold', rotation=90, labelpad=20, color='black')
        self.Ax.axhline(y=threshold, color='grey', linestyle='--', alpha=0.5, linewidth=2)
        threshold_handle = Line2D([0], [0], color='grey', linestyle='--', alpha=1, linewidth=2, label=f'$\\mathbf{{Thresh:\\ {threshold}}}$')
        legend_handles.append(threshold_handle)

        #match cluster to color and add to legend
        for cluster in sorted(cluster_to_color.keys()):
            color = cluster_to_color[cluster]
            legend_handles.append(mpatches.Patch(color=color, label=f"$\\mathbf{{Cluster\\ {cluster}}}$ :\n{cluster_counts[cluster]} datasets"))

        #Add legend
        self.Ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.02, 1),
                    borderaxespad=0, fontsize="small", handleheight=3, handlelength=3,
                    title_fontsize="medium")

        self.dendroPlot.tight_layout()

        if len(sorted_clusters) >= 10:
            self.dendroTitle.setText(f"<html><span style='color: red; font-weight: bold; font-size: 14px;'>Warning: Colors may not match clusters due to large number of clusters, please check cluster content tab</span><br>\
                                    <span style='font-weight: bold; font-size: 14px;'>Dendrogram with threshold: <span style='font-weight: bold;'>{threshold}</span>, the largest cluster is cluster \
                                    <span style='color: blue;'>{sorted_clusters[0][0]} with {sorted_clusters[0][1]} datasets</span>\
                                    The second largest cluster is cluster <span style='color: blue; font-weight: bold;'>{sorted_clusters[1][0]} with {sorted_clusters[1][1]} datasets</span></html>")
        elif len(sorted_clusters) >= 2 and len(sorted_clusters) < 10:
            self.dendroTitle.setText(f"<html><span style='color: green; font-weight: bold; font-size: 14px;'>Found less than 10 clusters, colors should match</span><br>\
                                    <span style='font-weight: bold; font-size: 14px;'>Dendrogram with threshold: <span style='color: blue;'>{threshold}</span>, the largest cluster is cluster <span style='color: blue;'>{sorted_clusters[0][0]} with {sorted_clusters[0][1]} datasets</span>\
                                    The second largest cluster is cluster <span style='color: blue;'>{sorted_clusters[1][0]} with {sorted_clusters[1][1]} datasets</span></html>")
        else:
            if len(sorted_clusters) == 1 and sorted_clusters[0][1] > 1:
                self.dendroTitle.setText(f"<html><span style='color: black; font-weight: bold;'>Dendrogram with threshold: {threshold}, the largest cluster is cluster {sorted_clusters[0][0]} with {sorted_clusters[0][1]} datasets</span></html>")
            else:
                self.dendroTitle.setText(f"<html><span style='color: black; font-weight: bold;'>Dendrogram with threshold: {threshold}, no clusters with more than 1 datasets found</span></html>")

        self.dendroCanvas.draw()


    #function to show the contents of merged cluster in a new tab, it will get some data from the Plotdrogram function, so it should be called after PlotDendrogram
    def showMergedCluster(self, ccClusterfile:str=None, threshold:float=None, SelectedClusterStr:str=None):
        #Check if ccClusterfile and threshold are valid
        if ccClusterfile is None:
            self.clusterContentText.setText("No ccCluster file in Work dir")
            return None
                
        if threshold is None:
            self.clusterContentText.setText("No threshold provided")
            return None

        #Get FlatC and etiquets from the ccClusterfile
        _, Tree, etiquets, _ = self.setupCC_method(ccClusterfile)
        FlatC = hierarchy.fcluster(Tree, threshold, criterion='distance')
        counter=collections.Counter(FlatC)

        #create a list of selected clusters, need to have SelectedClusterStr.strip() != "", as it may pass a empty string, which will cause error when split
        if SelectedClusterStr is not None and SelectedClusterStr.strip() != "":
            #change the SelectedClusterStr to a list of integers, if it fails, use the largest cluster instead
            #the selected cluster number needs to be a int, or it cannot compare with cluster number (np.unique(FlatC) returns a list of int)
            try:
                SelectedCluster = [int(x) for x in SelectedClusterStr.split()]
            except ValueError:
                print(f"Warning: Invalid cluster numbers provided: {SelectedClusterStr}. Using the largest cluster instead.")
                SelectedCluster = [max(counter.items(), key=operator.itemgetter(1))[0]]
        else:
            print(f"No selected clusters provided, using the largest cluster instead.")
            SelectedCluster = [max(counter.items(), key=operator.itemgetter(1))[0]]

        #Check whether FlatC indices matches etiques
        CheckIndices = checkIndices(FlatC, etiquets)
        if not CheckIndices:
            self.clusterContentText.setText(f"<html><span style='color: red; font-weight: bold; font-size: 18px;'>Warning: Mismatch between dataset number in ccClusterlogtxt and FlatC, \
                                            please check the ccCluster file</span></html>")
        else:
            #Get the indices of the datasets in each cluster and sort them by size
            clusterAndDataset = {}
            for cluster in np.unique(FlatC):
                indice = np.where(FlatC == cluster)[0].tolist()
                #use list comprehension to create a list of tuples (index, path) for each dataset in the cluster, and check whether the index matches the dataset number in etiquets
                tmpList = [(i, etiquets[i][1]) if int(etiquets[i][0]) == i else print(f"Warning: Mismatch in cluster assignment for dataset {etiquets[i][0]} and {i}") for i in indice]
                clusterAndDataset[cluster] = tmpList
            sorted_clusterAndDataset = dict(sorted(clusterAndDataset.items(), key=lambda item: len(item[1]), reverse=True))

            #for debug use, commented out to avoid cluttering the output
            #print(f"SortedClusterToIndices: {sorted_clusterAndDataset}")
            #print(f"etiquets: {etiquets}")
            """for cluster, datasets in sorted_clusterAndDataset.items():
                print(f"{cluster} has {len(datasets)} datasets:")
                for data in datasets:
                    print(f"{data}")"""

            #set title for the cluster content text with HTML format. For non-input text should be OK, otherwise use plain text to avoid the input text being interpreted as HTML
            self.clustercontentTitle.setText(f"<html><span style='color: black; font-weight: bold; font-size: 14px;'>Cluster content with threshold: \
                                                        {self.threshold_widget.text().strip() if self.threshold_widget else ''}; \
                                                      The merged cluster and the corresponding datasets will be shown in <span style='color: blue;'>BLUE</span> color</span></html>")

            #add the cluster content to the text area with the selected clusters highlighted in blue
            content = ""
            for cluster, datasets in sorted_clusterAndDataset.items():
                if SelectedCluster is not None and cluster in SelectedCluster:
                    content += f"<span style='color: blue; font-weight: bold;'>Cluster {cluster} ({len(datasets)} datasets):</span><br>"
                    for data in datasets:
                        content += f"<span style='color: blue;'>data number: {data[0]}; data path: {data[1]}</span><br>"
                else:
                    content += f"<span style='font-weight: bold;'>Cluster {cluster} ({len(datasets)} datasets):</span><br>"
                    for data in datasets:
                        content += f"data number: {data[0]}; data path: {data[1]}<br>"
            self.clusterContentText.setHtml(content)

            #add the path of the datasets in the selected clusters to the clusterPathText area for copying
            path_content = ""
            for cluster, datasets in sorted_clusterAndDataset.items():
                if SelectedCluster is not None and cluster in SelectedCluster:
                    path_content += f"Cluster {cluster} ({len(datasets)} datasets):\n"
                    for data in datasets:
                        path_content += f"{data[1]}\n"
                else:
                    path_content += f"Cluster {cluster} ({len(datasets)} datasets):\n"
                    for data in datasets:
                        path_content += f"{data[1]}\n"
            self.clusterPathText.setPlainText(path_content)
            
