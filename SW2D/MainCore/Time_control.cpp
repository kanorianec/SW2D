#include "Raschet.h"
#include <time.h>

time_t mkgmtime(struct tm* tm)
{
#if defined(_WIN32) || defined(_WIN64)
    return _mkgmtime(tm);
#elif defined(__linux__)
    return timegm(tm);
#endif
}
      

void Raschet::SetStartTime(int Year, int Month, int Day, int Hour, int Minute, int Second)
{
	struct tm tmp = { 0 };
	tmp.tm_year = Year - 1900;
	tmp.tm_mon = Month - 1;
	tmp.tm_mday = Day;
	tmp.tm_hour = Hour;
	tmp.tm_min = Minute;
	tmp.tm_sec = Second;
	tmp.tm_isdst = 0;
	
	RaschetTime = mkgmtime(&tmp);
}
