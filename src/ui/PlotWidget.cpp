#include "ui/PlotWidget.h"

#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <QStringList>
#include <QToolTip>

#include <algorithm>
#include <cmath>

PlotWidget::PlotWidget(QWidget* parent) : QWidget(parent) {
    setMinimumSize(480, 340);
    setMouseTracking(true);
    setAttribute(Qt::WA_OpaquePaintEvent);
}

void PlotWidget::setSeries(std::vector<Series> series) {
    series_ = std::move(series);
    updateBounds();
    update();
}

void PlotWidget::setAxisLabels(QString xLabel, QString yLabel) {
    xLabel_ = std::move(xLabel);
    yLabel_ = std::move(yLabel);
    update();
}

void PlotWidget::setEmptyMessage(QString message) {
    emptyMessage_ = std::move(message);
    update();
}

QRectF PlotWidget::plotRect() const {
    return QRectF(74.0, 28.0, std::max(10, width() - 104), std::max(10, height() - 90));
}

QPointF PlotWidget::mapPoint(const QPointF& point, const QRectF& rect) const {
    const double x = rect.left() + (point.x() - xMin_) / (xMax_ - xMin_) * rect.width();
    const double y = rect.bottom() - (point.y() - yMin_) / (yMax_ - yMin_) * rect.height();
    return {x, y};
}

void PlotWidget::updateBounds() {
    if (series_.empty()) {
        xMin_ = -1.0; xMax_ = 1.0; yMin_ = 0.0; yMax_ = 1.0;
        return;
    }
    bool found = false;
    for (const auto& series : series_) {
        for (const auto& point : series.points) {
            if (!std::isfinite(point.x()) || !std::isfinite(point.y())) continue;
            if (!found) {
                xMin_ = xMax_ = point.x();
                yMin_ = yMax_ = point.y();
                found = true;
            } else {
                xMin_ = std::min(xMin_, point.x());
                xMax_ = std::max(xMax_, point.x());
                yMin_ = std::min(yMin_, point.y());
                yMax_ = std::max(yMax_, point.y());
            }
        }
    }
    if (!found) {
        xMin_ = -1.0; xMax_ = 1.0; yMin_ = 0.0; yMax_ = 1.0;
        return;
    }
    if (std::abs(xMax_ - xMin_) < 1e-12) {
        xMin_ -= 1.0;
        xMax_ += 1.0;
    }
    if (std::abs(yMax_ - yMin_) < 1e-12) {
        yMin_ -= 0.5;
        yMax_ += 0.5;
    }
    const double xPadding = (xMax_ - xMin_) * 0.04;
    const double yPadding = (yMax_ - yMin_) * 0.10;
    xMin_ -= xPadding;
    xMax_ += xPadding;
    yMin_ = std::min(0.0, yMin_ - yPadding * 0.25);
    yMax_ += yPadding;
}

QString PlotWidget::formatNumber(double value) const {
    if (std::abs(value) >= 10000.0 || (std::abs(value) > 0.0 && std::abs(value) < 0.001))
        return QString::number(value, 'e', 1);
    return QString::number(value, 'g', 4);
}

void PlotWidget::paintEvent(QPaintEvent*) {
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.fillRect(rect(), QColor("#fbfcff"));

    const QRectF area = plotRect();
    if (series_.empty()) {
        painter.setPen(QColor("#8790a5"));
        QFont font = painter.font();
        font.setPointSize(11);
        painter.setFont(font);
        painter.drawText(rect(), Qt::AlignCenter, emptyMessage_);
        return;
    }

    QFont labelFont = painter.font();
    labelFont.setPointSize(8);
    painter.setFont(labelFont);
    painter.setPen(QColor("#e6e9f2"));
    constexpr int tickCount = 6;
    for (int i = 0; i <= tickCount; ++i) {
        const double ratio = static_cast<double>(i) / tickCount;
        const double x = area.left() + ratio * area.width();
        const double y = area.bottom() - ratio * area.height();
        painter.drawLine(QPointF(x, area.top()), QPointF(x, area.bottom()));
        painter.drawLine(QPointF(area.left(), y), QPointF(area.right(), y));

        painter.setPen(QColor("#737b91"));
        const double xValue = xMin_ + ratio * (xMax_ - xMin_);
        const double yValue = yMin_ + ratio * (yMax_ - yMin_);
        painter.drawText(QRectF(x - 36, area.bottom() + 8, 72, 18),
                         Qt::AlignHCenter | Qt::AlignTop, formatNumber(xValue));
        painter.drawText(QRectF(2, y - 9, 62, 18),
                         Qt::AlignRight | Qt::AlignVCenter, formatNumber(yValue));
        painter.setPen(QColor("#e6e9f2"));
    }

    painter.setPen(QPen(QColor("#aeb5c5"), 1));
    painter.drawLine(area.bottomLeft(), area.bottomRight());
    painter.drawLine(area.topLeft(), area.bottomLeft());
    painter.setPen(QColor("#4b5368"));
    painter.drawText(QRectF(area.left(), area.bottom() + 32, area.width(), 18),
                     Qt::AlignCenter, xLabel_);
    painter.save();
    painter.translate(18, area.center().y());
    painter.rotate(-90);
    painter.drawText(QRectF(-area.height() / 2, -10, area.height(), 18),
                     Qt::AlignCenter, yLabel_);
    painter.restore();

    painter.save();
    painter.setClipRect(area.adjusted(-1, -1, 1, 1));
    for (const auto& series : series_) {
        if (series.points.empty()) continue;
        if (series.style == Style::Bars) {
            const double fallbackWidth = area.width() / std::max<std::size_t>(series.points.size(), 1) * 0.78;
            for (std::size_t i = 0; i < series.points.size(); ++i) {
                const QPointF mapped = mapPoint(series.points[i], area);
                double barWidth = fallbackWidth;
                if (series.points.size() > 1) {
                    const std::size_t neighbor = i + 1 < series.points.size() ? i + 1 : i - 1;
                    barWidth = std::abs(mapPoint(series.points[neighbor], area).x() - mapped.x()) * 0.72;
                }
                const double baseY = mapPoint({series.points[i].x(), 0.0}, area).y();
                QColor fill = series.color;
                fill.setAlpha(145);
                painter.setBrush(fill);
                painter.setPen(QPen(series.color, 1));
                painter.drawRoundedRect(QRectF(mapped.x() - barWidth / 2.0, mapped.y(),
                                               barWidth, baseY - mapped.y()), 2, 2);
            }
            continue;
        }

        QPainterPath path;
        path.moveTo(mapPoint(series.points.front(), area));
        for (std::size_t i = 1; i < series.points.size(); ++i) {
            const QPointF previous = mapPoint(series.points[i - 1], area);
            const QPointF current = mapPoint(series.points[i], area);
            if (series.style == Style::Steps) {
                path.lineTo(current.x(), previous.y());
                path.lineTo(current);
            } else {
                path.lineTo(current);
            }
        }
        if (series.filled) {
            QPainterPath fillPath = path;
            fillPath.lineTo(mapPoint({series.points.back().x(), 0.0}, area));
            fillPath.lineTo(mapPoint({series.points.front().x(), 0.0}, area));
            fillPath.closeSubpath();
            QColor fill = series.color;
            fill.setAlpha(30);
            painter.fillPath(fillPath, fill);
        }
        painter.setPen(QPen(series.color, 2.2));
        painter.setBrush(Qt::NoBrush);
        painter.drawPath(path);
    }
    painter.restore();

    double legendX = area.right() - 8;
    for (auto it = series_.rbegin(); it != series_.rend(); ++it) {
        const int textWidth = painter.fontMetrics().horizontalAdvance(it->name);
        legendX -= textWidth + 34;
        painter.setPen(QPen(it->color, 3));
        painter.drawLine(QPointF(legendX, area.top() + 7), QPointF(legendX + 18, area.top() + 7));
        painter.setPen(QColor("#343b50"));
        painter.drawText(QPointF(legendX + 23, area.top() + 11), it->name);
        legendX -= 16;
    }

    if (mouseInside_ && area.contains(mousePosition_)) {
        painter.setPen(QPen(QColor(82, 92, 120, 80), 1, Qt::DashLine));
        painter.drawLine(QPointF(mousePosition_.x(), area.top()), QPointF(mousePosition_.x(), area.bottom()));
        const double x = xMin_ + (mousePosition_.x() - area.left()) / area.width() * (xMax_ - xMin_);
        QStringList values;
        values << QStringLiteral("x = %1").arg(formatNumber(x));
        for (const auto& series : series_) {
            if (series.points.empty()) continue;
            const auto nearest = std::min_element(series.points.begin(), series.points.end(),
                [x](const QPointF& a, const QPointF& b) {
                    return std::abs(a.x() - x) < std::abs(b.x() - x);
                });
            values << QStringLiteral("%1 = %2").arg(series.name, formatNumber(nearest->y()));
        }
        QToolTip::showText(mapToGlobal(mousePosition_ + QPoint(12, 12)), values.join('\n'), this);
    }
}

void PlotWidget::mouseMoveEvent(QMouseEvent* event) {
    mousePosition_ = event->position().toPoint();
    mouseInside_ = true;
    update();
}

void PlotWidget::leaveEvent(QEvent*) {
    mouseInside_ = false;
    QToolTip::hideText();
    update();
}
