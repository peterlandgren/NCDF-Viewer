NCDFViewer:
	python setup.py py2app
	cp -R /opt/local/lib/Resources/qt_menu.nib dist/NCDFViewer.app/Contents/Resources/
	cp -R qt.conf dist/NCDFViewer.app/Contents/Resources/
