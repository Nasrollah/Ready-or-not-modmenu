//////////////////////////////////
// Generated header: TypeTraits.h
// Forwards to the appropriate code
// that works on the detected compiler
// Generated on Mon Sep 30 23:14:48 2002
//////////////////////////////////

#ifdef LOKI_USE_REFERENCE
#	include "../../loki/TypeTraits.h"
#else
#	if (__INTEL_COMPILER)
#		include "../../loki/TypeTraits.h"
#	elif (__MWERKS__)
#		include "../../loki/TypeTraits.h"
#	elif (__BORLANDC__ >= 0x560)
#		include "../Borland/TypeTraits.h"
#	elif (_MSC_VER >= 1301)
#		include "../../loki/TypeTraits.h"
#	elif (_MSC_VER >= 1300)
#		include "../MSVC/1300/TypeTraits.h"
#	elif (_MSC_VER >= 1200)
#		include "../MSVC/1200/TypeTraits.h"
#	else
#		include "../../loki/TypeTraits.h"
#	endif
#endif
