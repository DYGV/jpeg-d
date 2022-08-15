import std.stdio : writef, writefln, writeln, File;
import jpeg_d;

void main() {
    string file_name = "sample.jpeg";
    file_name = "image.jpg";
    YCbCr img = jpeg_d.decode_jpeg(file_name);
    print_component(img);
}

void print_component(YCbCr img) {
    for(int i_c=0; i_c<img.component.length; i_c++){
        ubyte* component = &img.component[i_c].component[0];
        const int w = img.component[i_c].width;
        const int h = img.component[i_c].height;
        "w = ".writeln(w);
        "h = ".writeln(h);
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                writef("%3d ", component[w * i + j]);
            }
            writeln;
        }
    }
}

