#pragma once

#include <igl/colormap.h>
#include <imgui.h>

void drawColormap(igl::ColorMapType colormap){
    // drawing code : https://github.com/ocornut/imgui/issues/3606
    // widget: https://gist.github.com/ocornut/51367cc7dfd2c41d607bb0acfa6caf66
    ImVec2 size(ImGui::GetContentRegionAvailWidth(), 20.0f);
    ImGui::InvisibleButton("canvas", size);
    ImVec2 p0 = ImGui::GetItemRectMin();
    ImVec2 p1 = ImGui::GetItemRectMax();
    ImDrawList* draw_list = ImGui::GetWindowDrawList();
    draw_list->PushClipRect(p0, p1);

    float sx = 1.f / 16.f;
    float sy = 1.f / 5.f;
    Eigen::VectorXd gradient(100);
    for (int i=0; i<100; i++) gradient(i) = (double) i / 100.0;
    Eigen::MatrixXd rgb;
    igl::colormap(colormap, gradient, true, rgb);

    for (int i=0; i<100; i++)
    {
        float tx = (float) i / 100.0;
        ImVec2 c((tx + 0.5f * sx), (0 + 0.5f * sy));
        float k = 0.5f;
        draw_list->AddRectFilled(
            ImVec2(p0.x + (c.x - k * sx) * size.x, p0.y),
            ImVec2(p0.x + (c.x + k * sx) * size.x, p0.y + size.y),
            IM_COL32(255*rgb(i,0), 255*rgb(i,1), 255*rgb(i,2), 255));

        if (i == 100/2){
            draw_list->AddRectFilled(
                ImVec2(p0.x + (c.x - k * sx) * size.x, p0.y),
                ImVec2(p0.x + (c.x + k * sx) * size.x, p0.y + size.y),
                IM_COL32(0, 0, 0, 255));
        }
    }
    draw_list->PopClipRect();    
}