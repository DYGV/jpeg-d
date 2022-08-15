module jpeg_d.decode;

import std.stdio : writef, writefln, writeln, File;
import std.file : read;
import std.array : array;
import std.algorithm : sum;

import std.bitmanip : bitfields;
import std.math : cos, PI, round;
import std.conv : to;
import std.range : repeat;

// dfmt off
int[8][8] zig_zag_order = [
    [ 0,  1,  5,  6, 14, 15, 27, 28],
    [ 2,  4,  7, 13, 16, 26, 29, 42],
    [ 3,  8, 12, 17, 25, 30, 41, 43],
    [ 9, 11, 18, 24, 31, 40, 44, 53],
    [10, 19, 23, 32, 39, 45, 52, 54],
    [20, 22, 33, 38, 46, 51, 55, 60],
    [21, 34, 37, 47, 50, 56, 59, 61],
    [35, 36, 48, 49, 57, 58, 62, 63],
];
// dfmt on

align(1) struct APP {
    align(1) ubyte[5] JFIF = [0x4A, 0x46, 0x49, 0x46, 0x0];
    ushort _version;
    ubyte unit;
    ushort horizontal_resolution;
    ushort vertical_resolution;
    ubyte horizontal_thumbs_size;
    ubyte vertical_thumbs_size;
    ubyte[] img_thumbs;
}

align(1) struct SOS_C {
align(1):
    ubyte identify;
    mixin(bitfields!(ubyte, "ac_table_num", 4, ubyte, "dc_table_num", 4));
}

align(1) struct SOS {
align(1):
    ubyte num_of_component;
    SOS_C[3] component;
    ubyte q_start_num;
    ubyte q_end_num;
    mixin(bitfields!(ubyte, "prev_shift_quantity", 4, ubyte, "curr_shift_quantity",
            4));
}

align(1) struct DHT {
align(1):
    mixin(bitfields!(ubyte, "identify", 4, ubyte, "code", 4));
    ubyte[16] dist_num;
    ubyte[256] huff_val;

    ushort[] huff_code;
    ubyte[16] val_ptr;
    int[16] min_code;
    int[16] max_code;
}

align(1) struct DQT {
align(1):
    mixin(bitfields!(ubyte, "precision", 4, ubyte, "identify", 4));
    ubyte[64] table;
}

align(1) struct SOF_C {
align(1):
    ubyte identify;
    mixin(bitfields!(ubyte, "horizontal_sampling_factor", 4, ubyte,
            "vertical_sampling_factor", 4));
    ubyte q_table_selector;
}

align(1) struct SOF0 {
align(1):
    ubyte sample_precision;
    ushort image_row;
    ushort image_col;
    ubyte num_of_component;
    SOF_C[3] component;

    int h_max;
    int v_max;
}

struct INFO_C {
    int dc_pred;
    ubyte[] component;
    int width;
    int height;
}

struct YCbCr {
    INFO_C[3] component;
}

YCbCr g_info;
File file;

void copy_segment_data(T)(T* _struct, ubyte* data, int size) {
    ubyte* buf = cast(ubyte*)_struct;
    for (int i = 0; i < size; i++) {
        buf[i] = data[i];
    }
}

void parse_APP(APP* app, ubyte* data, int size) {
    copy_segment_data(app, data, size);
}

void parse_DHT(DHT* dht, ubyte* data, int size) {
    int table_offset = 0;
    copy_segment_data(dht, data, size);

    while (table_offset < size - 4) {
        table_offset += dht.dist_num.to!(int[]).sum + 17;
        ubyte[] huffsize;
        for (int i = 0; i < dht.dist_num.length; i++) {
            huffsize ~= i.repeat(dht.dist_num[i]).array.to!(ubyte[]);
        }
        int code = 0;
        int get_bit = 1;
        // huff size分コードを求める
        for (int i = 0; i < huffsize.length; i++) {
            // bit数が増加したらcodeをシフトする
            if (huffsize[i] > get_bit) {
                code = code << (huffsize[i] - get_bit);
            }
            // 記録用
            get_bit = huffsize[i];
            // huffsizeビット分取り出す
            dht.huff_code ~= cast(ushort)(code & (int.max >> huffsize[i]));
            code++;
        }
        int j = 0;
        for (int i = 0; i < 16; i++) {
            if (dht.dist_num[i] == 0) {
                dht.max_code[i] = -1;
            } else {
                dht.val_ptr[i] = cast(ubyte)j;
                dht.min_code[i] = dht.huff_code[j];
                j = j + dht.dist_num[i] - 1;
                dht.max_code[i] = dht.huff_code[j];
                j += 1;
            }
        }
    }
}

void swap_lower_upper_byte(ushort* data) {
    *data = ((*data >> 8) | (*data & 0x00ff));
}

void parse_SOF0(SOF0* sof, ubyte* data, int size) {
    copy_segment_data(sof, data, size);
    swap_lower_upper_byte(&sof.image_col);
    swap_lower_upper_byte(&sof.image_row);
}

void parse_DQT(DQT* dqt, ubyte* data, int size) {
    copy_segment_data(dqt, data, size);
}

void parse_SOS(SOS* sos, ubyte* data, int size) {
    copy_segment_data(sos, data, size);
}

DHT getDHT(DHT[] dht, int index, bool is_dc) {
    for (int i = 0; i < dht.length; i++) {
        if (is_dc && dht[i].identify == index && dht[i].code == 0) {
            return dht[i];
        }
        if (!is_dc && dht[i].identify == index && dht[i].code == 1) {
            return dht[i];
        }
    }
    return dht[0];
}

void calc_h_v_max(ref SOF0 sofn) {
    for (int i = 0; i < sofn.num_of_component; i++) {
        if (sofn.h_max < sofn.component[i].horizontal_sampling_factor) {
            sofn.h_max = sofn.component[i].horizontal_sampling_factor;
        }
        if (sofn.v_max < sofn.component[i].vertical_sampling_factor) {
            sofn.v_max = sofn.component[i].vertical_sampling_factor;
        }
    }
}

YCbCr decode(SOF0* sofn, SOS sos, DHT[] dht, DQT[] dqt) {
    const int xf = sofn.h_max << 3;
    const int yf = sofn.v_max << 3;
    const int mcux = (sofn.image_col + xf - 1) / xf;
    const int mcuy = (sofn.image_row + yf - 1) / yf;
    const int x = mcux * xf;
    const int y = mcuy * yf;

    for (int i = 0; i < sofn.num_of_component; i++) {
        const int w = mcux * sofn.component[i].horizontal_sampling_factor * 8;
        const int h = mcuy * sofn.component[i].vertical_sampling_factor * 8;
        g_info.component[i].component = new ubyte[w * h];
        g_info.component[i].width = w;
        g_info.component[i].height = h;
    }

    for (int row = 0; row < mcuy; row++) {
        for (int col = 0; col < mcux; col++) {
            decode_mcu(col, row, *sofn, sos, dht, dqt);
        }
    }
    return g_info;
}

void print_component(int index) {
    ubyte* component = &g_info.component[index].component[0];
    const int w = g_info.component[index].width;
    const int h = g_info.component[index].height;
    "w = ".writeln(w);
    "h = ".writeln(h);
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            writef("%3d ", component[w * i + j]);
        }
        writeln;
    }
}

void decode_mcu(int mcux, int mcuy, SOF0 sofn, SOS sos, DHT[] dht, DQT[] dqt) {
    for (int i = 0; i < sofn.num_of_component; i++) {
        const int Td = sos.component[i].dc_table_num;
        const int Ta = sos.component[i].ac_table_num;
        const int Tq = sofn.component[i].q_table_selector;
        DHT dht_dc = get_dht(Td, dht, true);
        DHT dht_ac = get_dht(Ta, dht, false);
        DQT _dqt = get_dqt(Tq, dqt);
        decode_mcu_component(mcux, mcuy, i,
                sofn.component[i].horizontal_sampling_factor,
                sofn.component[i].vertical_sampling_factor, dht_dc, dht_ac, _dqt);
    }
}

DHT get_dht(int index, DHT[] dht, bool is_dc) {
    for (int i = 0; i < dht.length; i++) {
        if (is_dc) {
            if ((dht[i].code == 0) && (dht[i].identify == index)) {
                return dht[i];
            }
        } else {
            if ((dht[i].code == 1) && (dht[i].identify == index)) {
                return dht[i];
            }
        }
    }
    return dht[0];
}

DQT get_dqt(int index, DQT[] dqt) {
    for (int i = 0; i < dqt.length; i++) {
        if (dqt[i].identify == index) {
            return dqt[i];
        }
    }
    return dqt[0];
}

void decode_mcu_component(int mcux, int mcuy, int index, int H, int V,
        DHT dc, DHT ac, DQT dqt) {
    int[8][8] F;
    int[8][8] f;

    for (int row = 0; row < V; row++) {
        for (int col = 0; col < H; col++) {
            int[64] ZZ;
            decode_entropy_coded_data_unit(index, dc, ac, &ZZ[0]);
            dequantize(dqt, &ZZ[0]);
            rearrange_zig_zag_order(&ZZ[0], F);
            IDCT(F, f);
            level_shift(f);
            copy_data_unit(mcux * H + col, mcuy * V + row, index, f);
        }
    }
}

void copy_data_unit(int dux, int duy, int index, int[8][8] f) {
    INFO_C* _component = &g_info.component[0];
    const int stride = _component[index].width;
    ubyte* component = &_component[index].component[0];
    const int offset = duy * 8 * stride + dux * 8;

    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            ubyte* p = component + offset + stride * i;
            p[j] = f[i][j].to!ubyte;
        }
    }
}

void IDCT(int[8][8] F, ref int[8][8] f) {
    double N = 8.0;

    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++) {
            double val = 0;
            for (int v = 0; v < 8; v++) {
                for (int u = 0; u < 8; u++) {
                    double cu, cv;
                    if (u == 0) {
                        cu = 0.707106781;
                    } else {
                        cu = 1.0;
                    }
                    if (v == 0) {
                        cv = 0.707106781;
                    } else {
                        cv = 1.0;
                    }
                    val += cu * cv * F[v][u] * cos(
                            (2.0 * x + 1.0) * u * PI / (2.0 * N)) * cos(
                            (2.0 * y + 1.0) * v * PI / (2.0 * N));
                }
            }
            f[y][x] = round((2.0 / N) * val).to!int;
        }
    }
}

void level_shift(ref int[8][8] f) {
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            f[i][j] += 128;
            if (f[i][j] < 0) {
                f[i][j] = 0;
            }
            if (f[i][j] > 255) {
                f[i][j] = 255;
            }
        }
    }
}

void rearrange_zig_zag_order(int* ZZ, ref int[8][8] F) {
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            F[i][j] = ZZ[zig_zag_order[i][j]];
        }
    }
}

void dequantize(DQT dqt, int* ZZ) {
    for (int i = 0; i < 64; i++) {
        ZZ[i] *= dqt.table[i];
    }
}

ubyte nextbyte() {
    ubyte[1] val;
    get_bytes(val);
    return val[0];
}

int nextbit() {
    static int cnt = 0;
    static ubyte b = 0;
    int bit = 0;

    if (cnt == 0) {
        b = nextbyte();
        cnt = 8;

        if (b == 0xFF) {
            nextbyte();
        }
    }

    bit = b >> 7;
    cnt = cnt - 1;
    b = cast(ubyte)(b << 1);

    return bit;
}

int decode(DHT dht, int* code, int* size) {
    int _code = 0;
    int i = 0;
    for (i = 0; i < 16; i++) {
        _code = (_code << 1) + nextbit();
        if (_code <= dht.max_code[i]) {
            break;
        }
    }
    int j = dht.val_ptr[i];
    j = j + _code - dht.min_code[i];
    return dht.huff_val[j];
}

int decode_dc_coefficient(DHT dht) {
    int code, size;
    const int T = decode(dht, &code, &size);
    const int DIFF = receive(T);
    return extend(DIFF, T);
}

int receive(int ssss) {
    int v = 0;

    for (int i = 0; i < ssss; i++) {
        v = (v << 1) + nextbit();
    }
    return v;
}

int extend(int V, int T) {
    int Vt = (1 << (T - 1));

    if (V < Vt) {
        Vt = (-1 << T) + 1;
        V = V + Vt;
    }

    return V;
}

void decode_entropy_coded_data_unit(int index, DHT dc, DHT ac, int* ZZ) {
    ZZ[0] = decode_dc_coefficient(dc) + g_info.component[index].dc_pred;
    g_info.component[index].dc_pred = ZZ[0];
    decode_ac_coefficients(ac, &ZZ[1]);
}

void decode_ac_coefficients(DHT dht, int* ZZ) {
    int K = 0;
    int code = 0;
    int size = 0;

    while (1) {
        const int RS = decode(dht, &code, &size);
        int SSSS = RS % 16;
        const int RRRR = RS >> 4;
        const int R = RRRR;

        if (SSSS == 0) {
            if (R == 15) {
                K += 16;
                if (K >= 62) {
                    break;
                }
            } else {
                break;
            }
        } else {
            K += R;
            ZZ[K] = receive(SSSS);
            ZZ[K] = extend(ZZ[K], SSSS);
            if (K >= 62) {
                break;
            } else {
                K += 1;
            }
        }
    }
}

void get_bytes(ubyte[] buffer) {
    file.rawRead(buffer);
}

YCbCr decode_jpeg(string file_name) {
    file = File(file_name, "rb");
    ubyte[2] marker;
    get_bytes(marker);
    if (marker[1] != 0xd8) {
        throw new Exception("error SOI");
    }
    APP[] app;
    DQT[] dqt;
    DHT[] dht;
    SOF0 sof0;
    SOS sos;
    while (1) {
        get_bytes(marker);
        if (marker[0] != 0xff) {
            throw new Exception("error: marker");
        }
        // EOI
        if (marker[1] == 0xd9) {
            break;
        }
        ubyte m = marker[1];
        ubyte[2] seg_size;
        get_bytes(seg_size);
        int size = ((seg_size[0] << 8) | seg_size[1]) - 2;
        ubyte[] seg_data = new ubyte[size];
        get_bytes(seg_data);

        switch (m) {
        case 0xdb: // 量子化テーブル
            DQT* dqt_ = new DQT();
            parse_DQT(dqt_, seg_data.ptr, size);
            dqt ~= *dqt_;
            break;

        case 0xc0: // フレームタイプ0(ベースライン)
            parse_SOF0(&sof0, seg_data.ptr, size);
            break;

        case 0xc4: // ハフマン符号のテーブル
            DHT* dht_ = new DHT();
            parse_DHT(dht_, seg_data.ptr, size);
            dht ~= *dht_;
            break;

        case 0xda: // スキャン
            parse_SOS(&sos, seg_data.ptr, size);
            calc_h_v_max(sof0);
            decode(&sof0, sos, dht, dqt);
            break;

        case 0xe0:
            APP* app_ = new APP();
            parse_APP(app_, seg_data.ptr, size);
            app ~= *app_;
            break;

        default:
            "ignored marker as not supported (0xFF%X)".writefln(m);
            break;
        }
    }
    return g_info;
}

void print_dqt(DQT[] dqt) {
    for (int i = 0; i < dqt.length; i++) {
        "DQT".writeln;
        "precision=%d, identify=%d".writefln(dqt[i].precision, dqt[i].identify);
        for (int j = 0; j < 8; j++) {
            for (int k = 0; k < 8; k++) {
                "%d, ".writef(dqt[i].table[(j * 8) + k]);
            }
            writeln;
        }
        "----------------------------".writeln;
    }
}

void print_dht(DHT[] dht) {
    for (int i = 0; i < dht.length; i++) {
        "DHT".writeln;
        "code class: ".writeln(dht[i].code);
        "identify  : ".writeln(dht[i].identify);
        "length:  code=%d,val=%d".writefln(dht[i].huff_code.length,
                dht[i].huff_val.length);
        "HUFFCODE: ".writeln(dht[i].huff_code);
        "HUFFVAL: ".writeln(dht[i].huff_val);
        "VALPTR: ".writeln(dht[i].val_ptr);
        "MAXCODE: ".writeln(dht[i].max_code);
        "MINCODE: ".writeln(dht[i].min_code);

        "----------------------------".writeln;
    }
}

void print_sof0(SOF0 sof0) {
    "SOF0".writeln;
    "row=%d, col=%d\nnumber_of_component=%d".writefln(sof0.image_row,
            sof0.image_col, sof0.num_of_component);
    for (int i = 0; i < sof0.num_of_component; i++) {
        "component_identify=%d: horizontal_sampling_factor=%d, vertical_sampling_factor=%d, DQT_identify=%d"
            .writefln(sof0.component[i].identify,
                    sof0.component[i].horizontal_sampling_factor,
                    sof0.component[i].vertical_sampling_factor,
                    sof0.component[i].q_table_selector);
    }
}

void print_sos(SOS sos) {
    "SOS".writeln;
    "num_of_component: ".writeln(sos.num_of_component);
    for (int i = 0; i < sos.num_of_component; i++) {
        "identify=%d, dc_huffman_table_num=%d, ac_huffman_table_num=%d"
            .writefln(sos.component[i].identify,
                    sos.component[i].dc_table_num, sos.component[i].ac_table_num);
    }
    writeln;

    "quantize start num: %d\nquantize end num: %d\nprev shift quantity%d\ncurr shift quantity:%d"
        .writefln(sos.q_start_num, sos.q_end_num,
                sos.prev_shift_quantity, sos.curr_shift_quantity);
}

void print_file_info(DQT[] dqt, DHT[] dht, SOF0 sof0, SOS sos) {
    print_dqt(dqt);
    "----------------------------".writeln;
    print_dht(dht);
    "----------------------------".writeln;
    print_sof0(sof0);
    "----------------------------".writeln;
    print_sos(sos);
    "----------------------------".writeln;
}
