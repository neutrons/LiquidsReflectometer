# How to batch reduce several runs

A batch reduction script is installed in the instrument shared folder. On an analysis computer (or lrac.sns.gov), start a terminal and do the following:

```
cd /SNS/REF_L/shared
```

```
python3 batch_reduce.py <IPTS> <first run> <last run>
```

Where you replace `<IPTS>` with your full IPTS number, `<first run>` and `<last run>` with the first and last run numbers of the range of runs you want to reduce.

For example, your command may look like this:
```
python3 batch_reduce.py IPTS-20406 178195 178216
```
The script will then reduce each run and post the results on https://monitor.sns.gov


# How to process a set of direct beam to produce scaling factors

A scaling factor script is installed in the instrument shared folder. On an analysis computer (or lrac.sns.gov), start a terminal and do the following:

```
cd /SNS/REF_L/shared
```

```
python3 process_sf.py <incident medium> <first run> <last run> <cfg file name>
```

Where you replace `<incident medium>` is the name of the incident medium you want to use, `<first run>` and `<last run>` with the first and last run numbers of the range of runs you want to use,
and `<cfg file name>>` is the full file path of your new scaling factor file.

For example, your command may look like this:
```
python process_sf.py Si 178195 178216 /SNS/REF_L/shared/autoreduce/sf_178195_Si2InDiam.cfg
```
