from goes2go.data import goes_timerange

from datetime import datetime
import pandas as pd

## Dates may be specified as datetime, pandas datetimes, or string dates
## that pandas can interpret.

## Specify start/end time as a panda-parsable string
start = '2024-07-19 00:00'
end   = '2024-07-23 00:00'

g = goes_timerange(start, end,
                   satellite='goes18',
                   product='ABI-L1b-Rad',
                   return_as='filelist')


