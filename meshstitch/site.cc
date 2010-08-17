#include "site.h"

Site::Site(void)
:n(0),orient(-1),l1(NULL),l2(NULL),th(NULL),ev(NULL)
{
};

Site::~Site(void)
{
  delete[] ev;
};

Site::Site (const Site& rhs)
:n(rhs.n),orient(rhs.orient),l1(rhs.l1),l2(rhs.l2),th(rhs.th),ev(rhs.ev)
{
}

