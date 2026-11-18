load("@rules_cc//cc:defs.bzl", "cc_library")
load("@rules_foreign_cc//foreign_cc:defs.bzl", "cmake")

cc_library(
    name = "glm",
    hdrs = glob([
        "glm/**/*.hpp",
        "glm/**/*.inl",
        "glm/**/*.h",
    ], allow_empty = True),
    includes = ["."],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "imgui",
    srcs = [
        "imgui.cpp",
        "imgui_draw.cpp",
        "imgui_tables.cpp",
        "imgui_widgets.cpp",
    ],
    hdrs = [
        "imgui.h",
        "imconfig.h",
        "imgui_internal.h",
        "imstb_rectpack.h",
        "imstb_textedit.h",
        "imstb_truetype.h",
    ],
    includes = ["."],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "tracy",
    srcs = ["public/TracyClient.cpp"],
    hdrs = glob([
        "public/**/*.h",
        "public/**/*.hpp",
    ], allow_empty = True),
    textual_hdrs = glob([
        "public/**/*.cpp",
    ], allow_empty = True),
    includes = ["public"],
    defines = ["TRACY_ENABLE"],
    linkopts = [
        "-pthread",
        "-ldl",
    ],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "imgui_sfml",
    srcs = ["imgui-SFML.cpp"],
    hdrs = glob(["*.h"], allow_empty = True),
    includes = ["."],
    deps = [
        "@imgui_src//:imgui",
        "@sfml_src//:sfml",
    ],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "json",
    hdrs = ["single_include/nlohmann/json.hpp"],
    strip_include_prefix = "single_include/",
    visibility = ["//visibility:public"],
)

    

filegroup(
    name = "all_srcs",
    srcs = glob(["**"]),
    visibility = ["//visibility:public"],
)
        
cmake(
    name = "sfml",
    lib_source = "@sfml_src//:all_srcs",
    out_static_libs = [
        "libsfml-graphics-s.a",
        "libsfml-window-s.a",
        "libsfml-system-s.a",
        "libsfml-audio-s.a",
        "libsfml-network-s.a",
    ],
    visibility = ["//visibility:public"],
)
    
    
