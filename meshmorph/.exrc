if &cp | set nocp | endif
map ' `
vmap [% [%m'gv``
vmap ]% ]%m'gv``
vmap a% [%v]%
let s:cpo_save=&cpo
set cpo&vim
nmap gx <Plug>NetrwBrowseX
nnoremap <silent> <Plug>NetrwBrowseX :call netrw#NetBrowseX(expand("<cWORD>"),0)
let &cpo=s:cpo_save
unlet s:cpo_save
set autoindent
set background=dark
set backspace=eol,start
set backupdir=/home/jed/.vim-backup,.
set cedit=<C-F>
set cinoptions=g0t0+3sc2s(0*80
set cmdheight=2
set cmdwinheight=10
set directory=/home/jed/.vim-temp,.
set display=lastline,uhex
set expandtab
set exrc
set fileencodings=ucs-bom,utf-8,default,latin1
set helplang=en
set hlsearch
set ignorecase
set incsearch
set laststatus=2
set lazyredraw
set listchars=tab:>-,eol:$,trail:+
set maxmem=16384
set maxmemtot=16384
set printoptions=paper:letter
set ruler
set runtimepath=~/.vim,/var/lib/vim/addons,/usr/share/vim/addons,/usr/share/vim/vimfiles,/usr/share/vim/vim70,/usr/share/vim/vimfiles/after,/usr/share/vim/addons/after,/var/lib/vim/addons/after,~/.vim/after
set scrolloff=2
set secure
set shiftwidth=2
set shortmess=aOtI
set showbreak=+++\ \ \ \ \ 
set showcmd
set showmatch
set sidescrolloff=3
set smartcase
set smarttab
set softtabstop=2
set nostartofline
set suffixes=.bak,~,.swp,.o,.info,.aux,.log,.dvi,.bbl,.blg,.brf,.cb,.ind,.idx,.ilg,.inx,.out,.toc
set tildeop
set viminfo='50,\"1000,h
set visualbell
set whichwrap=b,s,[,]
set wildignore=*.o,*.class,*.a
set wildmode=longest,list
set window=50
set winminheight=0
set winminwidth=20
" vim: set ft=vim :
