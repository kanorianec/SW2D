#include "Raschet.h"
#include <time.h>

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
	
	RaschetTime = _mkgmtime(&tmp);
}
