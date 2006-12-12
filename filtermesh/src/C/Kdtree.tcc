template <typename T>
Kdtree<T>::Kdtree(int pnd, int pmaxlevel)
	: nd(pnd),
	  maxlevel(pmaxlevel),
	  fsize(0)
{
	assertx(nd>0 && maxlevel>0);
	root=new KdtreeNode();	// creates illegal tree because axis==-1
	if (getenv("KDFSIZE")) fsize=atof(getenv("KDFSIZE"));
}

template <typename T>
Kdtree<T>::~Kdtree()
{
	if (KDTREE_STATS) {
		STAT(SKDdepth);
		recDepth(root,SKDdepth,0);
	}
	recClear(root);
}

template <typename T>
void Kdtree<T>::clear()
{
	if (KDTREE_STATS) { STAT(SKDdepth); recDepth(root,SKDdepth,0); }
	recClear(root);
	root=new KdtreeNode();
}

template <typename T>
void Kdtree<T>::enter(T id, const double bb[][2])
{
	double avgel=0;
	if (fsize) {
		ForIndex(i,nd) { avgel+=bb[i][1]-bb[i][0]; } EndFor;
		avgel/=nd;
	}
	SArray<double> aval(nd); ForIndex(i,nd) { aval[i]=.5; } EndFor;
	recEnter(id,bb,root,aval,0,.5,0,avgel);
}

template <typename T>
int Kdtree<T>::search(double bb[][2], CBFUNC cbfunc, KdtreeNode *n) const
{
	if (!n) n=root;
	int nelemvis=0;
	int ret=n->axis<0 ? 0 : recSearch(n,n,bb,cbfunc,nelemvis); // check empty tree
	SKDsearchnel+=nelemvis;
	return ret;
}

template <typename T>
void Kdtree<T>::recClear(KdtreeNode* n)
{
	if (!n) return;
	recClear(n->l);
	recClear(n->h);
	while (!n->stack.empty()) Kdentry::specialdelete(n->stack.pop(),nd);
	delete n;
}

template <typename T>
void Kdtree<T>::recDepth(KdtreeNode* n, Stat& stat, int depth) const
{
	if (!n) return;
	ForStack(n->stack,Kdentry*,e) { (void)e; stat+=depth; } EndFor;
	recDepth(n->l,stat,depth+1);
	recDepth(n->h,stat,depth+1);
}

template <typename T>
void Kdtree<T>::recEnter(T id, const double bb[][2],
		      KdtreeNode* n, SArray<double>& aval,
		      int level, double inc, int axis, double avgel)
{
	for (;;) {
		const double val=aval[axis];
		if (n->axis<0) n->axis=axis,n->val=val;
		if (!axis) {
			if (++level==maxlevel) break;
			inc*=.5;
		}
		int wl=bb[axis][0]<=val; // want l
		int wh=bb[axis][1]>=val; // want h
		if (wl && wh) {		 // single recursion
			if (!fsize || avgel>=inc*fsize)
				break; // small enough
			SArray<double> naval(nd);
			ForIndex(i,nd) { naval[i]=aval[i]; } EndFor;
			naval[axis]-=inc;
			if (!n->l) n->l=new KdtreeNode();
			recEnter(id,bb,n->l,naval,level,inc,(axis+1)%nd,avgel);
			if (!n->h) n->h=new KdtreeNode();
			n=n->h; aval[axis]+=inc;
		} else if (wl) { // go down l
			if (!n->l) n->l=new KdtreeNode();
			n=n->l; aval[axis]-=inc;
		} else if (wh) { // go down h
			if (!n->h) n->h=new KdtreeNode();
			n=n->h; aval[axis]+=inc;
		} else {
			assertnever("");
		}
		axis=(axis+1)%nd;
	}
	Kdentry* e=Kdentry::specialnew(nd); e->id = Conv<T>::e(id);
	ForIndex(i,nd) { e->bb[i][0]=bb[i][0]; e->bb[i][1]=bb[i][1]; } EndFor;
	n->stack.push(e);
}

// lca==lowest common ancestor
template <typename T>
int Kdtree<T>::recSearch(KdtreeNode* n, KdtreeNode* lca,
			 double bb[][2], CBFUNC cbfunc, int &nelemvis) const
{
        int i;
	for (;;) {
		ForStack(n->stack,Kdentry*,e) {
			nelemvis++;
			// cast for DECCXX
			const double (*ebb)[2]=e->bb;
//			for (int i=0;i<nd;i++)
			for (i=0;i<nd;i++)
				if (ebb[i][0]>=bb[i][1] || ebb[i][1]<=bb[i][0])
					break;
			if (i<nd) continue;
			if (cbfunc(Conv<T>::d(e->id),bb,lca)==2)
				return 1;
		} EndFor;
		const int axis=n->axis;
		const double val=n->val;
		int wl=n->l && bb[axis][0]<val;
		int wh=n->h && bb[axis][1]>val;
		if (wl && wh) { // single recursion
			if (recSearch(n->h,lca,bb,cbfunc, nelemvis)) return 1;
			// bb may have changed, test again
			if (!(bb[axis][0]<val)) return 0;
			n=n->l;
		} else if (wl) { // no recursion
			if (lca==n) lca=n->l;
			n=n->l;
		} else if (wh) { // no recursion
			if (lca==n) lca=n->h;
			n=n->h;
		} else {
			return 0;
		}
	}
}

template <typename T>
void Kdtree<T>::recPrint(KdtreeNode* n, int l) const
{
        int i;
//	for (int i=0;i<l;i++) std::cerr << " ";
	for (i=0;i<l;i++) std::cerr << " ";
	if (!n) { std::cerr << "<nil>\n"; return; }
	SHOWF("partition of axis %d along %g <<\n",n->axis,n->val);
	ForStack(n->stack,Kdentry*,e) {
		for (i=0;i<l;i++) std::cerr << " ";
		std::cerr << Conv<int>::d(e->id) << "\n";
	} EndFor;
	for (i=0;i<l;i++) std::cerr << " ";
	std::cerr << ">>\n";
	recPrint(n->l,l+1);
	recPrint(n->h,l+1);
}
