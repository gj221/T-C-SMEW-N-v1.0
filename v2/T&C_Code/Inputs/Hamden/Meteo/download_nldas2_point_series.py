"""NLDAS-2 primary forcing (NLDAS_FORA0125_H v2.0) at the Hamden parcel as
hourly point time series, one CSV per variable, through NASA Harmony's
Giovanni time-series service.  Credentials from ~/.netrc."""
import netrc, datetime as dt, os, shutil, sys
from harmony import Client, Collection, Request, WKT
from harmony.config import Environment

OUT = r"C:\NLDAS\Hamden\giovanni"
POINT = "POINT(-72.866 41.437)"
START, STOP = dt.datetime(2014, 12, 31, 0), dt.datetime(2026, 1, 1, 23, 59, 59)
SUFFIX = ""
if len(sys.argv) == 4:   # start, stop (yyyy-mm-ddTHH), file suffix
    START, STOP, SUFFIX = dt.datetime.fromisoformat(sys.argv[1]), dt.datetime.fromisoformat(sys.argv[2]), "_" + sys.argv[3]
VARS = ["Tair", "Qair", "PSurf", "Wind_E", "Wind_N", "SWdown", "Rainf"]

login, _, password = netrc.netrc().authenticators("urs.earthdata.nasa.gov")
c = Client(auth=(login, password), env=Environment.PROD)
os.makedirs(OUT, exist_ok=True)
for v in VARS:
    target = os.path.join(OUT, f"NLDAS_FORA0125_H_2.0_{v}{SUFFIX}.csv")
    if os.path.exists(target):
        print(v, "already downloaded", flush=True)
        continue
    req = Request(collection=Collection(id="C2033151148-GES_DISC"), spatial=WKT(POINT),
                  temporal={"start": START, "stop": STOP},
                  variables=[f"NLDAS_FORA0125_H_2_0_{v}"], format="text/csv")
    job = c.submit(req)
    print(v, "job", job, flush=True)
    c.wait_for_processing(job, show_progress=False)
    tmp = os.path.join(OUT, "_tmp_" + v)
    os.makedirs(tmp, exist_ok=True)
    files = [f.result() for f in c.download_all(job, directory=tmp, overwrite=True)]
    if len(files) != 1:
        sys.exit(f"{v}: expected one result file, got {files}")
    shutil.move(files[0], target)
    shutil.rmtree(tmp)
    print(v, "->", target, os.path.getsize(target), "bytes", flush=True)
print("ALL DONE", flush=True)
