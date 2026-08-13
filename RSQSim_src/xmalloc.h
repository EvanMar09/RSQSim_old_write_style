#ifndef XMALLOC_H
#define XMALLOC_H

#include <stdlib.h>
void *xmalloc(size_t size);
void *xcalloc(size_t nmemb, size_t size);
void *xrealloc(void *ptr, size_t size);
void xfree(void *ptr);

#endif /* ndef XMALLOC_H */
