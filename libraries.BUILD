load("@rules_cc//cc:defs.bzl", "cc_library")
load("@rules_foreign_cc//foreign_cc:defs.bzl", "cmake")

config_setting(
    name = "tracy_enabled",
    define_values = {"tracy": "false"},
)


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
    name = "imgui_sfml",
    srcs = ["imgui-SFML.cpp"],
    hdrs = glob(["*.h"], allow_empty = True),
    includes = ["."],
    deps = [
        "@//:sfml",
        "@imgui_src//:imgui",
    ],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "json",
    hdrs = ["single_include/nlohmann/json.hpp"],
    strip_include_prefix = "single_include/",
    visibility = ["//visibility:public"],
)

cc_library(
    name = "tracy",
    srcs = select({
        ":tracy_enabled": ["public/TracyClient.cpp"],
        "//conditions:default": [],
    }),
    hdrs = glob([
        "public/**/*.h",
        "public/**/*.hpp",
    ], allow_empty = True),
    textual_hdrs = glob([
        "public/**/*.cpp",
    ], allow_empty = True),
    includes = ["public"],
    defines = select({
        ":tracy_enabled": ["TRACY_ENABLE"],
        "//conditions:default": [],
    }),
    linkopts = select({
        ":tracy_enabled": ["-pthread", "-ldl"],
        "//conditions:default": [],
    }),
    visibility = ["//visibility:public"],
)
    
    
