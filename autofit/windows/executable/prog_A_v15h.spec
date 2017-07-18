# -*- mode: python -*-

block_cipher = None

mkl_dlls = [('C:\\Python27\\Scripts\\mk2_avx.dll', ''),
 ('C:\\Python27\\Scripts\\mk2_def.dll', ''),
 ('C:\\Python27\\Scripts\\mk2_vml_avx.dll', ''),
 ('C:\\Python27\\Scripts\\mk2_vml_def.dll', '')]

a = Analysis(['prog_A_v15h.py'],
             pathex=['C:\\triples\\v15 - Copy'],
             binaries=mkl_dlls,
             datas=[('C:\\Python27\\tcl\\tcl8.5', 'lib\\tcl8.5'), ('C:\\Python27\\tcl\\tk8.5', 'lib\\tk8.5')],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='prog_A_v15h',
          debug=False,
          strip=False,
          upx=True,
          console=True )
