from __future__ import print_function, absolute_import

__author__ = "Rita Giordano, Gianluca Santoni, Tiankun Zhou"
__copyright__ = "Copyright 2015-2026"
__credits__ = ["Rita Giordano, Gianluca Santoni, Tiankun Zhou, Alexander Popov"]
__license__ = ""
__version__ = "2.0"
__maintainer__ = "Tiankun Zhou"
__email__ = "tiankun.zhou@esrf.fr"
__status__ = "Beta"


from PySide6 import QtWidgets
import os, sys

from .textSummary import generateLogSummary
#a class to generate the results widget.
#will be used as a tab in the main window

#not used anymore
class resultsSummary(QtWidgets.QWidget):
    def __init__(self):
        QtWidgets.QWidget.__init__(self)
        self.resultsLayout= QtWidgets.QVBoxLayout(self)
        self.setGeometry(100, 100, 500, 500)
        self.Workdir= os.getcwd()
#the title of the results summary
        self.Title=QtWidgets.QLabel(self)
        self.Title.setText('Summary of the results')
        self.resultSummary= QtWidgets.QTextEdit()
        self.resultsLayout.addWidget(self.Title)
        self.resultsLayout.addWidget(self.resultSummary)
        self.setText()

    def setText(self):
#the text edit to sum up results
        generateLogSummary(self.Workdir+'/.cc_cluster.log')
        text = open(self.Workdir+'/.cc_summary.txt').read()
        self.resultSummary.setPlainText(text)
        print(text)

def main():
#Can be launched alone to viasualize the results
#without loading all the program
    app = QtWidgets.QApplication(sys.argv)
    ex = resultsSummary()
    ex.show()
    sys.exit(app.exec())

if __name__== '__main__':
    main()
