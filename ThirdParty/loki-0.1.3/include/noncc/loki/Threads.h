///////////////////////////////
// Generated header: Threads.h
// Forwards to the appropriate code
// that works on the detected compiler
// Generated on Mon Sep 30 23:14:48 2002
///////////////////////////////

#ifdef LOKI_USE_REFERENCE
#	include "../../loki/Threads.h"
#else
#	if (__INTEL_COMPILER)
#		include "../../loki/Threads.h"
#	elif (__MWERKS__)
#		include "../../loki/Threads.h"
#	elif (__BORLANDC__ >= 0x560)
#		include "../Borland/Threads.h"
#	elif (_MSC_VER >= 1301)
#		include "../../loki/Threads.h"
#	elif (_MSC_VER >= 1300)
#		include "../MSVC/1300/Threads.h"
#	elif (_MSC_VER >= 1200)
#		include "../MSVC/1200/Threads.h"
#	else
#		include "../../loki/Threads.h"
#	endif
#endif
