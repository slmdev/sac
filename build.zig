const std = @import("std");
const fs = std.fs;
const mem = std.mem;
const path = std.fs.path;
const zcc = @import("compile_commands");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});
    const exe = b.addExecutable(.{
        .name = "sac",
        .target = target,
        .optimize = optimize,
    });

    collectCppFiles(b, exe) catch |err| {
        std.debug.print("File collection failed: {s}\n", .{@errorName(err)});
        @panic("Construction Aborted");
    };

    exe.linkLibCpp();
    b.installArtifact(exe);

    var targets = std.ArrayList(*std.Build.Step.Compile).init(b.allocator);
    targets.append(exe) catch @panic("OOM");
    zcc.createStep(b, "cdb", targets.toOwnedSlice() catch @panic("OOM"));
}

fn collectCppFiles(b: *std.Build, exe: *std.Build.Step.Compile) !void {
    const src_dir = "src";
    var cpp_files = std.ArrayList([]const u8).init(b.allocator);
    defer cpp_files.deinit();

    try collectFilesRecursive(
        b.allocator,
        src_dir,
        &cpp_files,
    );

    exe.addCSourceFiles(.{
        .files = cpp_files.items,
        .flags = &[_][]const u8{ "-std=c++23", "-fexperimental-library" },
    });
}

fn collectFilesRecursive(
    allocator: std.mem.Allocator,
    dir_path: []const u8,
    files: *std.ArrayList([]const u8),
) !void {
    var dir = try fs.cwd().openDir(dir_path, .{ .iterate = true });
    defer dir.close();

    var iter = dir.iterate();
    while (try iter.next()) |entry| {
        const full_path = try path.join(allocator, &[_][]const u8{ dir_path, entry.name });
        switch (entry.kind) {
            .directory => try collectFilesRecursive(allocator, full_path, files),
            .file => {
                if (mem.endsWith(u8, entry.name, ".cpp")) {
                    try files.append(full_path);
                }
            },
            else => {},
        }
    }
}
