"""Read x64 PE imports (including delay imports) without third-party modules."""
from pathlib import Path
import struct


def imports(path: Path):
    data = path.read_bytes()
    pe = struct.unpack_from('<I', data, 0x3c)[0]
    if data[:2] != b'MZ' or data[pe:pe+4] != b'PE\0\0':
        raise ValueError(f'Not a PE file: {path}')
    machine, nsections = struct.unpack_from('<HH', data, pe+4)
    optional_size = struct.unpack_from('<H', data, pe+20)[0]
    optional = pe+24
    if machine != 0x8664 or struct.unpack_from('<H', data, optional)[0] != 0x20b:
        raise ValueError(f'Expected x64 PE32+: {path}')
    image_base = struct.unpack_from('<Q', data, optional+24)[0]
    sections = []
    for i in range(nsections):
        size, rva, raw_size, offset = struct.unpack_from('<IIII', data, optional+optional_size+40*i+8)
        sections.append((rva, max(size, raw_size), offset))

    def file_offset(address):
        for rva, size, offset in sections:
            if rva <= address < rva+size:
                return offset+address-rva
        raise ValueError(f'Unmapped PE address in {path}: {address:x}')

    result = set()
    for index, width, name_offset in ((1, 20, 12), (13, 32, 4)):
        rva, size = struct.unpack_from('<II', data, optional+112+index*8)
        if not rva:
            continue
        offset = file_offset(rva)
        for pos in range(offset, offset+size, width):
            descriptor = data[pos:pos+width]
            if not any(descriptor):
                break
            address = struct.unpack_from('<I', descriptor, name_offset)[0]
            if index == 13 and not struct.unpack_from('<I', descriptor)[0] & 1:
                address -= image_base
            start = file_offset(address)
            result.add(data[start:data.index(b'\0', start)].decode('ascii').lower())
    return sorted(result)


def is_system(name):
    return name.startswith(('api-ms-win-', 'ext-ms-win-')) or name in {
        'kernel32.dll', 'kernelbase.dll', 'ntdll.dll', 'advapi32.dll', 'user32.dll',
        'gdi32.dll', 'shell32.dll', 'ole32.dll', 'oleaut32.dll', 'ws2_32.dll',
        'rpcrt4.dll', 'version.dll', 'winmm.dll', 'psapi.dll', 'bcrypt.dll',
        'secur32.dll', 'crypt32.dll', 'normaliz.dll', 'msvcrt.dll', 'ucrtbase.dll',
        'shlwapi.dll', 'comdlg32.dll', 'dbghelp.dll', 'iphlpapi.dll', 'imm32.dll',
        'imagehlp.dll', 'wintrust.dll',
    }
