set tabstop=4
set shiftwidth=4
set softtabstop=4

set expandtab
set smarttab

"Stops auto-line wrapping when typing.
set textwidth=0
set wrapmargin=0

set autoindent

set nu
set background=dark

set ignorecase
set smartcase
set incsearch

syntax enable

au BufRead,BufNewFile *.for let b:fortran_fixed_source=1
au BufRead,BufNewFile *.for let b:fortran_fixed_source=1

nmap <F1> :echo<CR>
imap <F1> <C-o>:echo<CR>

let fortran_have_tabs=1

