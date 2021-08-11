prefix := /home/controls/var/tmp/60Hz_Scans

check:
	diff standards/direct-39-25_Air_0att.csv $(prefix)/direct-39-25_Air_0att.csv
	diff standards/direct-39-25_Air.csv $(prefix)/direct-39-25_Air.csv
	diff standards/direct-39-25_Si2InDiam_0att.csv $(prefix)/direct-39-25_Si2InDiam_0att.csv
	diff standards/direct-39-25_Si2InDiam.csv $(prefix)/direct-39-25_Si2InDiam.csv
	diff standards/refl-39-25_0.2max.csv $(prefix)/refl-39-25_0.2max.csv
	diff standards/refl-39-25_0.2Qmax_inverted.csv $(prefix)/refl-39-25_0.2Qmax_inverted.csv
	diff standards/refl-39-25.csv $(prefix)/refl-39-25.csv
	diff standards/refl-39-25_inverted.csv $(prefix)/refl-39-25_inverted.csv

copy:
	cp standards/direct-39-25_Air_0att.csv $(prefix)/direct-39-25_Air_0att.csv
	cp standards/direct-39-25_Air.csv $(prefix)/direct-39-25_Air.csv
	cp standards/direct-39-25_Si2InDiam_0att.csv $(prefix)/direct-39-25_Si2InDiam_0att.csv
	cp standards/direct-39-25_Si2InDiam.csv $(prefix)/direct-39-25_Si2InDiam.csv
	cp standards/refl-39-25_0.2max.csv $(prefix)/refl-39-25_0.2max.csv
	cp standards/refl-39-25_0.2Qmax_inverted.csv $(prefix)/refl-39-25_0.2Qmax_inverted.csv
	cp standards/refl-39-25.csv $(prefix)/refl-39-25.csv
	cp standards/refl-39-25_inverted.csv $(prefix)/refl-39-25_inverted.csv

.PHONY: check copy
