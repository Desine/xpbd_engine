load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library")


cc_library(
    name = "sfml",
    linkopts = [
        "-lsfml-graphics",
        "-lsfml-window",
        "-lsfml-system",
        "-lsfml-audio",
        "-lsfml-network",
    ],
    visibility = ["//visibility:public"],
)

cc_binary(
    name = "exe",
    srcs = glob([
        "xpbd/**/*.cpp",
        "xpbd/**/*.h",
        "xpbd/**/*.hpp",
    ], allow_empty = True),
    deps = [
        ":sfml",
        "@glm_src//:glm",
        "@imgui_src//:imgui",
        "@imgui_sfml_src//:imgui_sfml",
        "@tracy_src//:tracy",
        "@nlohmann_json//:json",
    ],
    copts = [
        "-std=c++26",
        "-fopenmp",
    ],
    linkopts = [
        "-fopenmp",
        "-lX11",
        "-lXrandr", 
        "-lXi",
        "-lXcursor",
        "-lxcb",
        "-lxcb-image",
        "-ludev",
        "-lGL",
        "-lfreetype",
        "-lpthread",
        "-ldl",
    ],
)