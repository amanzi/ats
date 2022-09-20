import visit as v
import visit_rcParams as vrc
import datetime

def getTime(years=None, seconds=None, zero=None, time_format=None, round=None):
    assert not (years is None and seconds is None)
    
    if zero is None:
        zero = vrc.rcParams['time.zero']

    if time_format is None:
        time_format = vrc.rcParams['time.format']

    if years is not None:
        dt = datetime.timedelta(seconds=years*86400*365.25)
    else:
        dt = datetime.timedelta(seconds=seconds)

    time = zero + dt
    if round is not None:
        time = roundTime(time, round)

    return time.strftime(time_format)
    

def roundTime(time, round=None):
    """Rounds a time to the nearest Month (day 1)

    i.e. Jan 2, 10:00 am --> Jan 1
         Jan 29, 1:00 pm --> Feb 1    
    """
    if round is None:
        return time
    
    if round == "month":
        tprev = datetime.datetime(*(time.timetuple()[0:2]+(1,)))

        if time.timetuple()[1] == 12:
            yr_next = time.timetuple()[0]+1
            mo_next = 1
        else:
            yr_next = time.timetuple()[0]
            mo_next = time.timetuple()[1] + 1
        tnext = datetime.datetime(yr_next, mo_next, 1)

    assert tprev <= time <= tnext
    if time - tprev < tnext - time:
        tround = tprev
    else:
        tround = tnext
    print("Rounding: %s to %s"%(time.strftime("%c"), tround.strftime("%c")))
    return tround


def visitTime(round=None):
    v.Query("Time")
    return getTime(years=v.GetQueryOutputValue(), round=round)
